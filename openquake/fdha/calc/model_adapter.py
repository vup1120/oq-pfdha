# -*- coding: utf-8 -*-
"""
Legacy Model Adapter for FDHA calculations.

This module provides backward compatibility by wrapping existing FDHA models
to work with the new FDHAContext-based calculation pattern.
"""

import numpy as np
import inspect
import logging
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

# Near-field displacement floor: the smallest across-strike distance (km) fed to
# a distributed displacement regression that diverges as r -> 0 (Petersen 2011
# eq.18, ln r term). Fixed at half a 25-m Petersen cell (12.5 m). It exists only
# to tame the divergence, so it is a constant -- deliberately NOT the footprint
# z, which is a per-model occurrence-table cell size and can legitimately be
# hundreds of metres. Hard-coded and not exposed in job configuration.
# See docs/design/rupture_location_uncertainty.md (D7).
NEAR_FIELD_FLOOR_KM = 0.0125


class LegacyModelAdapter:
    """
    Wraps legacy FDHA models for use with FDHAContext.
    
    Provides backward compatibility during migration to context-based
    calculations. Automatically detects model type and handles parameter
    extraction from context objects.
    
    Example:
        adapter = LegacyModelAdapter(my_model, {'style': 'normal'})
        P_sr = adapter.compute_primary_sr(ctx, red_cfg)
    """
    
    def __init__(self, model: Any, model_params: Dict[str, Any] = None):
        """
        Initialize adapter for a model.
        
        Args:
            model: Model instance with get_prob() method
            model_params: Default parameters for model calls
        """
        self.model = model
        self.model_params = model_params or {}
        self._signature_cache: Dict[tuple, set] = {}
        
        # Auto-detect model type from class name
        class_name = model.__class__.__name__
        if 'PrimarySR' in class_name or 'PrimarySurfRup' in class_name:
            self.model_type = 'primary_sr'
        elif 'PrimaryFD' in class_name or 'PrimaryDispl' in class_name or 'PrimarySurfDispl' in class_name:
            self.model_type = 'primary_fd'
        elif 'SecondarySR' in class_name or 'SecondarySurfRup' in class_name:
            self.model_type = 'secondary_sr'
        elif 'SecondaryFD' in class_name or 'SecondaryDispl' in class_name or 'SecondarySurfDispl' in class_name:
            self.model_type = 'secondary_fd'
        else:
            self.model_type = 'unknown'
            logger.warning(f"Could not determine model type for: {class_name}")
    
    def _get_method_params(self, method_name: str) -> set:
        """
        Get valid parameter names for a method (cached).
        
        Args:
            method_name: Name of the method
            
        Returns:
            Set of valid parameter names
        """
        cache_key = (id(self.model.__class__), method_name)
        
        if cache_key not in self._signature_cache:
            method = getattr(self.model, method_name, None)
            if method is None:
                self._signature_cache[cache_key] = set()
            else:
                try:
                    sig = inspect.signature(method)
                    self._signature_cache[cache_key] = set(sig.parameters.keys())
                except (ValueError, TypeError):
                    self._signature_cache[cache_key] = set()
        
        return self._signature_cache[cache_key]
    
    def _call_safely(self, method_name: str, **kwargs) -> Any:
        """
        Call model method with only the parameters it accepts.
        
        Args:
            method_name: Name of method to call
            **kwargs: All possible parameters
            
        Returns:
            Method result or None on error
        """
        valid_params = self._get_method_params(method_name)
        filtered = {k: v for k, v in kwargs.items() if k in valid_params}
        
        try:
            return getattr(self.model, method_name)(**filtered)
        except Exception as e:
            logger.error(f"Error calling {self.model.__class__.__name__}.{method_name}: {e}")
            return None
    
    def _ctx_metrics(self, ctx: 'FDHAContext'):
        """(r, x_L, L) arrays for this model's declared reference-line method.

        Each FDHA model declares how multi-section (multiFaultSource)
        rupture distances must be measured via its MULTIFAULT_REFERENCE_LINE
        class attribute (ecs | lcp | segments); the context carries one
        metric set per method required by the configured models. For
        single-strand ruptures this is simply the canonical trace-based set.
        """
        method = getattr(self.model, 'MULTIFAULT_REFERENCE_LINE', 'lcp')
        return ctx.metrics_for(method)

    def compute_primary_sr(
        self,
        ctx: 'FDHAContext',
        red_cfg: Dict[str, Any]
    ) -> Optional[np.ndarray]:
        """
        Compute primary surface rupture probability.
        
        Primary SR is typically magnitude-dependent only, but some models
        may include epistemic uncertainty (MC samples).
        
        Args:
            ctx: FDHA context with rupture/site parameters
            red_cfg: MC reduction config {'method': 'median', 'q': 50}
            
        Returns:
            Array of shape (N,) or None on error
        """
        from openquake.fdha.calc.utils.probability import _reduce_mc
        
        N = len(ctx)
        
        # Style priority: model_params > derived from rake
        style = self.model_params.get('style')
        if style is None:
            style = ctx.style[0]  # Derived from rake angle
        
        # Handle models that don't support all styles (e.g., Youngs2003 only accepts 'all' or 'normal')
        model_name = self.model.__class__.__name__
        if 'Youngs2003' in model_name and style not in ('all', 'normal'):
            raise ValueError(
                f"Model {model_name} only supports style='all' or style='normal', "
                f"but got '{style}' (derived from rake={ctx.rake[0]}). "
                f"Please specify style explicitly in job.ini:\n"
                f"  [models.primary_surf_rup.parameters]\n"
                f"  style = all"
            )
        
        # Build kwargs; vs30 comes from the site collection (reference_vs30_value
        # when no per-site value is given) and is needed by e.g. Moss2013PrimarySR.
        kwargs = {
            'mag': float(ctx.mag[0]),
            'dip': float(ctx.dip[0]),
            'dip_mu': float(ctx.dip[0]),
            'seismothickness': self.model_params.get('seismothickness', 15.0),
            'vs30': float(ctx.vs30[0]),
            'style': style,
            **{k: v for k, v in self.model_params.items() if k != 'style'},
        }
        
        def _reduced_scalar(res):
            res_red = _reduce_mc(
                res,
                method=red_cfg.get('method', 'median'),
                q=red_cfg.get('q', 50)
            )
            return float(np.atleast_1d(res_red).flat[0])

        # vs30-dependent models (e.g. Moss 2013) give a different P_sr per
        # site when the site collection carries heterogeneous vs30 values;
        # evaluate once per unique vs30 instead of broadcasting site 0's.
        unique_vs30 = np.unique(ctx.vs30)
        if unique_vs30.size > 1 and 'vs30' in self._get_method_params('get_prob'):
            out = np.zeros(N, dtype=np.float64)
            for v in unique_vs30:
                result = self._call_safely('get_prob', **{**kwargs, 'vs30': float(v)})
                if result is None:
                    return None
                out[ctx.vs30 == v] = _reduced_scalar(result)
            return out

        # Call model
        result = self._call_safely('get_prob', **kwargs)

        if result is None:
            return None

        # Broadcast to all sites
        return np.full(N, _reduced_scalar(result), dtype=np.float64)
    
    def compute_primary_fd(
        self,
        ctx: 'FDHAContext',
        displacements: np.ndarray,
        red_cfg: Dict[str, Any]
    ) -> Optional[np.ndarray]:
        """
        Compute primary fault displacement probability.
        
        Uses fully vectorized operations - no per-displacement fallback loops.
        
        Args:
            ctx: FDHA context
            displacements: Target displacement levels (m)
            red_cfg: MC reduction config
            
        Returns:
            Array of shape (N, D) or None on error
        """
        from openquake.fdha.calc.utils.probability import _to_sites_x_displ
        
        N = len(ctx)
        D = len(displacements)
        
        # Style priority: model_params > derived from rake
        style = self.model_params.get('style')
        if style is None:
            style = ctx.style[0]  # Derived from rake angle
        
        # Validate style for models with limited support
        model_name = self.model.__class__.__name__
        if 'Youngs2003' in model_name and style not in ('all', 'normal'):
            raise ValueError(
                f"Model {model_name} only supports style='all' or style='normal', "
                f"but got '{style}' (derived from rake={ctx.rake[0]}). "
                f"Please specify style explicitly in job.ini:\n"
                f"  [models.primary_surf_displ.parameters]\n"
                f"  style = all"
            )
        
        # Build kwargs with vectorized arrays; x_L follows the model's
        # declared multi-fault reference line (e.g. Chiou2025 -> ECS).
        _r_sel, x_L_sel, _L_sel = self._ctx_metrics(ctx)
        kwargs = {
            'mag': float(ctx.mag[0]),
            'd': displacements,
            'X_L_ratio': x_L_sel,
            'x_L': x_L_sel,
            'style': style,
            **{k: v for k, v in self.model_params.items() if k != 'style'},
        }

        # Vectorized call - no fallback loops
        result = self._call_safely('get_prob', **kwargs)
        
        if result is None:
            logger.warning(f"compute_primary_fd returned None for {self.model.__class__.__name__}")
            return np.zeros((N, D), dtype=np.float64)
        
        arr = np.asarray(result)
        
        # Handle different output shapes
        if arr.shape == (N, D):
            return arr.astype(np.float64)
        elif arr.shape == (D, N):
            return arr.T.astype(np.float64)
        elif arr.ndim >= 2:
            return _to_sites_x_displ(arr, N, D, red_cfg)
        elif arr.ndim == 1:
            # 1D array - try to reshape
            if arr.size == N * D:
                return arr.reshape(N, D).astype(np.float64)
            elif arr.size == D:
                # Same probability for all sites - broadcast
                return np.broadcast_to(arr.reshape(1, D), (N, D)).copy().astype(np.float64)
            elif arr.size == N:
                # Per-site scalar - broadcast to all displacements
                return np.broadcast_to(arr.reshape(N, 1), (N, D)).copy().astype(np.float64)
        
        # Scalar or unknown shape - use _to_sites_x_displ for normalization
        return _to_sites_x_displ(arr, N, D, red_cfg)
    
    def compute_secondary_sr(
        self,
        ctx: 'FDHAContext',
        red_cfg: Dict[str, Any]
    ) -> np.ndarray:
        """
        Compute secondary (distributed) surface rupture probability.
        
        Uses fully vectorized operations - no per-site fallback loops.
        
        Args:
            ctx: FDHA context
            red_cfg: MC reduction config
            
        Returns:
            Array of shape (N,)
        """
        from openquake.fdha.calc.utils.probability import _reduce_mc
        
        N = len(ctx)
        
        # Style priority: model_params > derived from rake
        style = self.model_params.get('style')
        if style is None:
            style = ctx.style[0]  # Derived from rake angle
        
        # Validate style for models with limited support
        model_name = self.model.__class__.__name__
        if 'Youngs2003' in model_name and style not in ('all', 'normal'):
            raise ValueError(
                f"Model {model_name} only supports style='all' or style='normal', "
                f"but got '{style}' (derived from rake={ctx.rake[0]}). "
                f"Please specify style explicitly in job.ini:\n"
                f"  [models.secondary_surf_rup.parameters]\n"
                f"  style = all"
            )
        
        # Build kwargs with vectorized arrays; r follows the model's declared
        # multi-fault reference line (e.g. Visini2025 -> nearest segment).
        r_sel, _x_L_sel, _L_sel = self._ctx_metrics(ctx)
        kwargs = {
            'mag': float(ctx.mag[0]),
            'r': r_sel,
            'rx': ctx.rx,
            's': r_sel * 1000.0,
            'style': style,
            **{k: v for k, v in self.model_params.items() if k != 'style'},
        }

        # Ensure version is string if present (Youngs2003SecondarySR expects string)
        if 'version' in kwargs:
            kwargs['version'] = str(kwargs['version'])
        
        # Vectorized call - no fallback to per-site loops
        result = self._call_safely('get_prob', **kwargs)
        
        if result is None:
            logger.warning(f"compute_secondary_sr returned None for {self.model.__class__.__name__}")
            return np.zeros(N, dtype=np.float64)
        
        arr = np.asarray(result)
        
        # Handle different output shapes
        if arr.ndim == 0:
            # Scalar result - broadcast to all sites
            return np.full(N, float(arr), dtype=np.float64)
        
        if arr.ndim == 1:
            if arr.size == N:
                # Per-site results
                return arr.astype(np.float64)
            elif arr.size == 1:
                # Single value - broadcast
                return np.full(N, float(arr[0]), dtype=np.float64)
            else:
                # Reduce MC samples and broadcast
                reduced = _reduce_mc(
                    arr,
                    method=red_cfg.get('method', 'median'),
                    q=red_cfg.get('q', 50)
                )
                return np.full(N, float(np.atleast_1d(reduced).flat[0]), dtype=np.float64)
        
        if arr.ndim == 2:
            # (N, n_mc) or (n_mc, N) shape - reduce MC dimension
            if arr.shape[0] == N:
                reduced = _reduce_mc(
                    arr,
                    method=red_cfg.get('method', 'median'),
                    q=red_cfg.get('q', 50)
                )
                result_arr = np.atleast_1d(reduced)
                if result_arr.size == N:
                    return result_arr.astype(np.float64)
                else:
                    return np.full(N, float(result_arr.flat[0]), dtype=np.float64)
            elif arr.shape[1] == N:
                reduced = _reduce_mc(
                    arr.T,
                    method=red_cfg.get('method', 'median'),
                    q=red_cfg.get('q', 50)
                )
                result_arr = np.atleast_1d(reduced)
                if result_arr.size == N:
                    return result_arr.astype(np.float64)
                else:
                    return np.full(N, float(result_arr.flat[0]), dtype=np.float64)
            else:
                # Unknown shape - reduce and broadcast
                reduced = _reduce_mc(
                    arr,
                    method=red_cfg.get('method', 'median'),
                    q=red_cfg.get('q', 50)
                )
                return np.full(N, float(np.atleast_1d(reduced).flat[0]), dtype=np.float64)
        
        # Higher dimensional - reduce and broadcast
        reduced = _reduce_mc(
            arr,
            method=red_cfg.get('method', 'median'),
            q=red_cfg.get('q', 50)
        )
        return np.full(N, float(np.atleast_1d(reduced).flat[0]), dtype=np.float64)
    
    def compute_secondary_fd(
        self,
        ctx: 'FDHAContext',
        displacements: np.ndarray,
        red_cfg: Dict[str, Any],
    ) -> np.ndarray:
        """
        Compute secondary fault displacement probability.

        Args:
            ctx: FDHA context
            displacements: Target displacement levels (m)
            red_cfg: MC reduction config

        The near-field distance floor is the fixed NEAR_FIELD_FLOOR_KM
        constant; the distributed occurrence cell size lives with the
        secondary_surf_rup model's own pixel_size (FD logic tree). Neither
        is a caller-supplied parameter.

        Returns:
            Array of shape (N, D)
        """
        from openquake.fdha.calc.utils.probability import _to_sites_x_displ

        N = len(ctx)
        D = len(displacements)
        
        # Style priority: model_params > derived from rake
        style = self.model_params.get('style')
        if style is None:
            style = ctx.style[0]  # Derived from rake angle
        
        # Validate style for models with limited support
        model_name = self.model.__class__.__name__
        if 'Youngs2003' in model_name and style not in ('all', 'normal'):
            raise ValueError(
                f"Model {model_name} only supports style='all' or style='normal', "
                f"but got '{style}' (derived from rake={ctx.rake[0]}). "
                f"Please specify style explicitly in job.ini:\n"
                f"  [models.secondary_surf_displ.parameters]\n"
                f"  style = all"
            )
        
        # Build kwargs; r/x_L/L follow the model's declared multi-fault
        # reference line (e.g. Visini2025 -> nearest segment, raw GC2 x/L).
        r_sel, x_L_sel, L_sel = self._ctx_metrics(ctx)

        # Near-field floor (D7). Petersen (2011) eq.18 diverges as r -> 0; a
        # model declaring NEAR_FIELD_FLOOR == 'footprint_half' has the distance
        # fed to its displacement regression clamped to max(r, NEAR_FIELD_FLOOR_KM).
        # The floor is a FIXED 12.5 m (= half a 25-m Petersen cell), deliberately
        # decoupled from the footprint z: z can be a large model-specific cell
        # (e.g. 500 m selects a coarser occurrence table via pixel_size), and a
        # z/2 = 250 m floor would silently erase the near-trace distributed
        # hazard the model exists to produce -- it would push every on-trace
        # site out past the entire distributed zone. 12.5 m only tames ln(r),
        # nothing more. Hard-coded, not user-facing. The clamp lives here, at the
        # adapter boundary, so the model's get_prob stays paper-faithful. Bounded
        # models (Visini, Takao) declare nothing and are untouched.
        # See docs/design/rupture_location_uncertainty.md.
        if getattr(self.model, 'NEAR_FIELD_FLOOR', None) == 'footprint_half':
            r_sel = np.maximum(np.asarray(r_sel, dtype=np.float64),
                               NEAR_FIELD_FLOOR_KM)

        kwargs = {
            'mag': float(ctx.mag[0]),
            'd': displacements,
            'r': r_sel,
            'rx': ctx.rx,
            's': r_sel * 1000.0,
            'X_L_ratio': x_L_sel,
            'x_L': x_L_sel,
            'dip': ctx.dip,
            'L': L_sel,
            'style': style,
            **{k: v for k, v in self.model_params.items() if k != 'style'},
        }
        
        # Ensure percentile is string if present (Youngs2003SecondaryFD expects string)
        if 'percentile' in kwargs:
            kwargs['percentile'] = str(kwargs['percentile'])
        
        # Try vectorized call
        result = self._call_safely('get_prob', **kwargs)
        
        if result is not None:
            normalized = _to_sites_x_displ(result, N, D, red_cfg)
            return normalized
        
        # Return zeros on error
        logger.warning(f"compute_secondary_fd failed for {self.model.__class__.__name__}")
        return np.zeros((N, D), dtype=np.float64)