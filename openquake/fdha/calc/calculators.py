import os
import inspect
import logging
import numpy as np
from tqdm import tqdm
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.config_loader import load_config
from openquake.fdha.primary_surf_rup import *
from openquake.fdha.primary_surf_displ import *
from openquake.fdha.secondary_surf_rup import *
from openquake.fdha.secondary_surf_displ import *

logger = logging.getLogger(__name__)

class BaseFaultRuptureCalculator:
    """
    Base class for fault rupture probability calculations.
    """
    def __init__(self, config_path, source_model_paths, hdf5path='',
                 fault_sources=None, **converterparams):
        # Save config path (as absolute path) for path resolution
        import os
        if isinstance(config_path, (str, bytes, os.PathLike)):
            self.config_path = os.path.abspath(config_path)
        else:
            self.config_path = None
        # Load config
        self.config = self._load_configuration(config_path)
        # Parse source models via OQ-Engine, unless the caller pre-built them
        # (used by the logic-tree driver after applying NRML uncertainties).
        if fault_sources is None:
            self.fault_sources = parse_source_model_faults(
                source_model_paths, hdf5path, **converterparams,
            )
        else:
            self.fault_sources = dict(fault_sources)
        # Instantiate probability models
        self._initialize_models()
        # Initialize calculation parameters
        self._initialize_calculation_params()
        logger.debug(f"Initialized {self.__class__.__name__} with config {config_path}")

    @staticmethod
    def _load_configuration(config_path):
        if isinstance(config_path, (str, bytes, os.PathLike)):
            return load_config(config_path)
        else:
            raise TypeError(f"config_path must be a v5 canonical INI filepath, got {type(config_path)}")

    def _initialize_models(self):
        # Read model definitions from config
        models_cfg = self.config.get('models', {})
        self.primary_surf_rup_model = self._instantiate_model(models_cfg.get('primary_surf_rup'))
        self.primary_surf_displ_model = self._instantiate_model(models_cfg.get('primary_surf_displ'))
        self.secondary_surf_rup_model = self._instantiate_model(models_cfg.get('secondary_surf_rup'))
        self.secondary_surf_displ_model = self._instantiate_model(models_cfg.get('secondary_surf_displ'))
        # Multi-fault reference-line consistency: the principal FD model sets
        # the convention for the whole branch — the other models' distances
        # are measured against the same reference line so principal and
        # distributed hazard share one geometry (e.g. Chiou 2025 [ecs] pulls
        # a Petersen 2011 secondary onto the ECS line). The only exemption is
        # 'segments' (Visini 2025): its regression measured r to the nearest
        # rupturing section, which no smoothed line can represent.
        if self.primary_surf_displ_model is not None:
            _principal_line = getattr(self.primary_surf_displ_model,
                                      'MULTIFAULT_REFERENCE_LINE', 'lcp')
            for _m in (self.primary_surf_rup_model,
                       self.secondary_surf_rup_model,
                       self.secondary_surf_displ_model):
                if _m is not None and getattr(
                        _m, 'MULTIFAULT_REFERENCE_LINE', 'lcp') != 'segments':
                    _m.MULTIFAULT_REFERENCE_LINE = _principal_line
        logger.debug(f"Models: primary_surf_rup={self.primary_surf_rup_model}, primary_surf_displ={self.primary_surf_displ_model}, "
                     f"secondary_surf_rup={self.secondary_surf_rup_model}, secondary_surf_displ={self.secondary_surf_displ_model}")

    @staticmethod
    def _instantiate_model(model_cfg):
        if not model_cfg or 'type' not in model_cfg:
            return None
        model_type = model_cfg['type']
        params = model_cfg.get('parameters', {})
        
        # Look up class by name in current globals or imported modules.
        # A typo'd model name must fail the job loudly: returning None here
        # would silently zero the corresponding hazard contribution.
        try:
            model_class = globals()[model_type]
        except KeyError:
            raise ValueError(
                f"Unknown FDHA model class '{model_type}'; check the "
                f"uncertaintyModel / [models] configuration")

        # Check which parameters the model's __init__ accepts
        sig = inspect.signature(model_class.__init__)
        init_params = set(sig.parameters.keys()) - {'self'}

        # Only pass parameters that the constructor accepts
        constructor_params = {k: v for k, v in params.items() if k in init_params}

        # Constructor errors (e.g. an invalid scaling_model) must propagate:
        # swallowing them here used to convert a bad configuration into a
        # silently-zero hazard contribution.
        return model_class(**constructor_params)

    def call_model_safely(self, model, method, **kwargs):
        if model is None:
            logger.warning(f"No model for method {method}")
            return None
        try:
            func = getattr(model, method)
            sig = inspect.signature(func)
            # Pass only supported parameters
            filtered = {k: v for k, v in kwargs.items() if k in sig.parameters}
            return func(**filtered)
        except Exception as e:
            logger.error(f"Error calling {model}.{method}: {e}")
            return None

    def get_model_parameters(self, name):
        return self.config.get('models', {}).get(name, {}).get('parameters', {})

    @staticmethod
    def _validate_hazard_reduction(name, cfg):
        """Validate a Monte-Carlo reduction config for the hazard integral.

        Developer-facing; jobs should not set these keys. Only 'mean' and
        'median' are accepted here. The mean (default) is exact: expectation
        is linear, so reducing a model's epistemic/MC samples by their mean
        inside the rate sum reproduces the mean hazard curve (McGuire,
        Cornell & Toro, 2005). The median is a legacy central-estimate
        heuristic. 'percentile' is refused because a per-rupture quantile of
        exceedance probabilities is not a fractile of any hazard
        distribution — quantiles are only additive over the rate sum under
        comonotonicity (Dhaene et al., 2002). The percentile machinery in
        ``utils.probability`` is intentionally kept for non-integral uses
        (e.g. a future scenario calculator).
        """
        if isinstance(cfg, str):
            cfg = {'method': cfg}
        if not isinstance(cfg, dict):
            raise TypeError(
                f"[parameters] {name} must be a JSON object like "
                f'{{"method": "mean"}}, got {cfg!r}')
        method = str(cfg.get('method', 'mean')).lower()
        if method == 'percentile':
            raise ValueError(
                f'[parameters] {name} = {{"method": "percentile"}} is not '
                f'supported in hazard calculations: a quantile applied '
                f'inside the hazard integral does not produce a hazard '
                f'fractile of any kind. Use {{"method": "mean"}} — the '
                f'mean hazard curve, which incorporates within-model '
                f'epistemic uncertainty exactly. Fractiles of within-model '
                f'epistemic uncertainty are not currently supported; they '
                f'require propagating the model\'s posterior samples as '
                f'logic-tree realizations.')
        if method not in ('mean', 'median'):
            raise ValueError(
                f"[parameters] {name}: unknown reduction method "
                f"'{method}'; expected 'mean' or 'median'.")
        return {**cfg, 'method': method}

    def _initialize_calculation_params(self):
        """Initialize calculation parameters needed by calculate_fdha_hazard."""
        # Target displacements
        target_disp = self.config.get('parameters', {}).get('target_displacement', [0.001, 0.01, 0.1, 1.0, 10.0])
        self.target_displacements = np.array(target_disp, dtype=np.float64)
        
        # Reduction configs for models that return an internal MC/epistemic
        # sample dimension. The default is the mean, which is exact (the
        # expectation commutes with the hazard integral), so jobs never need
        # to set these keys; they are undocumented developer knobs.
        # secondary_sr_reduction inherits primary_sr_reduction when unset.
        self.p_sr_red_cfg = self._validate_hazard_reduction(
            'primary_sr_reduction',
            self.config.get('parameters', {}).get(
                'primary_sr_reduction', {'method': 'mean'}))
        self.s_sr_red_cfg = self._validate_hazard_reduction(
            'secondary_sr_reduction',
            self.config.get('parameters', {}).get(
                'secondary_sr_reduction', self.p_sr_red_cfg))
        
        # Principal/distributed split: r <= threshold -> primary (on-trace)
        # models, r > threshold -> secondary (distributed) models. Read from
        # [calculation] first, then [parameters]. Default 0.1 km.
        self.r_threshold_km = float(
            self.config.get('calculation', {}).get('r_threshold_km') or
            self.config.get('parameters', {}).get('r_threshold_km', 0.1)
        )

        # Near/far regime split used *inside* the secondary (Visini) SR Rank 2
        # Monte Carlo; distinct from r_threshold_km and not a model selector.
        # Default 0.2 km.
        self.near_far_threshold_km = float(
            self.config.get('parameters', {}).get('near_far_threshold_km', 0.2)
        )

        # Rupture-location uncertainty
        # (docs/design/rupture_location_uncertainty.md, section 3). W_p(r)
        # weights the principal contribution; G(r) weights the distributed
        # contribution per `combination_mode`. All read [calculation] first,
        # [parameters] as fallback, same idiom as r_threshold_km above.
        #
        #   r_sigma_km          two-sided mapping-accuracy sigma (Petersen
        #                       Tables 2-3). 0 = legacy boxcar; >0 = pinned
        #                       +/-n-sigma normal footprint mass.
        #   r_sigma_truncation  +/-n-sigma cut (Petersen p. 819 uses 2).
        #   site_footprint_m    footprint z (Petersen cell size); the sigma>0
        #                       window and the near-field-floor scale.
        self.r_sigma_km = float(
            self.config.get('calculation', {}).get('r_sigma_km') or
            self.config.get('parameters', {}).get('r_sigma_km', 0.0)
        )
        self.r_sigma_truncation = float(
            self.config.get('calculation', {}).get('r_sigma_truncation') or
            self.config.get('parameters', {}).get('r_sigma_truncation', 2.0)
        )
        self.site_footprint_m = float(
            self.config.get('calculation', {}).get('site_footprint_m') or
            self.config.get('parameters', {}).get('site_footprint_m', 25.0)
        )
        # combination_mode selects the distributed weight G(r):
        #   'additive'      G = 1 (Petersen et al. 2011 eq.1 + eq.2: principal
        #                   and distributed are independent and summed).
        #   'complementary' G = 1 - W_p (Youngs 2003 / Takao 2013 Fig.1 per-fault
        #                   either/or bookkeeping; the tool's historical split).
        # C2 flips this default from 'complementary' to 'additive'.
        _combination_mode = (
            self.config.get('calculation', {}).get('combination_mode') or
            self.config.get('parameters', {}).get('combination_mode') or
            'complementary'
        )
        self.combination_mode = str(_combination_mode).strip().lower()
        if self.combination_mode not in ('additive', 'complementary'):
            raise ValueError(
                f"[calculation] combination_mode: unknown value "
                f"'{self.combination_mode}'; expected 'additive' or "
                f"'complementary'.")

        # Case label for Visini models. The logic-tree branch typically sets
        # 'case' as a parameter of the secondary-model uncertaintyModel (it
        # travels with the Visini2025SecondarySR/FD branch, not [parameters]),
        # so fall back to those sections before defaulting.
        self.case_label = (
            self.config.get('parameters', {}).get('case')
            or self.get_model_parameters('secondary_surf_displ').get('case')
            or self.get_model_parameters('secondary_surf_rup').get('case')
            or 'case1'
        )

        # Union of the reference-line treatments the configured models
        # declare for multi-section (multiFaultSource) ruptures, via their
        # MULTIFAULT_REFERENCE_LINE class attribute — the FDHA analogue of
        # hazardlib collecting the union of the GMPEs' REQUIRES_DISTANCES.
        # The context maker computes one metric set per method in this union.
        _models = (self.primary_surf_rup_model, self.primary_surf_displ_model,
                   self.secondary_surf_rup_model, self.secondary_surf_displ_model)
        self.multifault_reference_lines = tuple(sorted(
            {getattr(m, 'MULTIFAULT_REFERENCE_LINE', 'lcp')
             for m in _models if m is not None} or {'lcp'}
        ))
        
        # Initialize model adapters
        from openquake.fdha.calc.model_adapter import LegacyModelAdapter
        self.adapters = {}
        
        if self.primary_surf_rup_model:
            self.adapters['primary_sr'] = LegacyModelAdapter(
                self.primary_surf_rup_model,
                self.get_model_parameters('primary_surf_rup')
            )
        
        if self.primary_surf_displ_model:
            self.adapters['primary_fd'] = LegacyModelAdapter(
                self.primary_surf_displ_model,
                self.get_model_parameters('primary_surf_displ')
            )
        
        if self.secondary_surf_rup_model:
            self.adapters['secondary_sr'] = LegacyModelAdapter(
                self.secondary_surf_rup_model,
                self.get_model_parameters('secondary_surf_rup')
            )
        
        if self.secondary_surf_displ_model:
            self.adapters['secondary_fd'] = LegacyModelAdapter(
                self.secondary_surf_displ_model,
                self.get_model_parameters('secondary_surf_displ')
            )
    
    def get_fdha_params(self):
        """
        Get the FDHA distance thresholds for the context maker.

        Returns a dict of two independent thresholds (km):
            - 'r_threshold_km': principal vs distributed split. Selects whether
              a site is handled by the primary (on-trace) or secondary
              (distributed) model family during hazard integration.
            - 'near_far_threshold_km': near vs far regime used *within* the
              secondary (Visini) SR Rank 2 along-strike Monte Carlo; it tunes
              the secondary computation rather than selecting the model family.
        """
        return {
            'r_threshold_km': self.r_threshold_km,
            'near_far_threshold_km': self.near_far_threshold_km,
            'multifault_reference_lines': self.multifault_reference_lines,
            # Rupture-location uncertainty (W_p / G); see
            # docs/design/rupture_location_uncertainty.md.
            'r_sigma_km': self.r_sigma_km,
            'r_sigma_truncation': self.r_sigma_truncation,
            'site_footprint_m': self.site_footprint_m,
            'combination_mode': self.combination_mode,
        }


class FaultRuptureProbabilityCalculator(BaseFaultRuptureCalculator):
    """
    Calculator for point-site fault rupture probabilities (hazard curves).
    """
    def __init__(self, config_path, source_model_paths, hdf5path='',
                 fault_sources=None, **converterparams):
        super().__init__(
            config_path, source_model_paths, hdf5path,
            fault_sources=fault_sources, **converterparams,
        )
        self.sitecol = None
        self._initialize_site()

    def _initialize_site(self):
        site_cfg = self.config.get('site_location', {})
        sites_list = site_cfg.get('sites_list')

        if sites_list:
            # Multi-site path. The actual hazard accumulation downstream
            # (FDHAContextMaker.get_ctx + _compute_rupture_contribution + np.add.at)
            # is numpy-vectorised across all N sites per rupture. The O(N) loop
            # here is a one-time initialisation cost, not per-rupture.
            n = len(sites_list)
            if n > 100:
                logger.warning(
                    f"Multi-site job has {n} sites (> 100). "
                    "Runtime scales with site count via the per-rupture distance loop "
                    "in VectorizedRuptureDistanceCalculator."
                )
            # SiteCollection([Site(...)]) is required (not from_points) so that
            # the 'vs30' field is present in the structured array — FDHAContextMaker
            # reads sitecol.vs30 and from_points omits that field.
            global_vs30 = site_cfg.get('vs30')
            oq_sites = []
            for s in sites_list:
                pt = Point(s['longitude'], s['latitude'])
                site_vs30 = s.get('vs30', global_vs30)
                kwargs = {'location': pt}
                if site_vs30 is not None:
                    kwargs['vs30'] = float(site_vs30)
                oq_sites.append(Site(**kwargs))
            self.sitecol = SiteCollection(oq_sites)
            logger.debug(f"Initialized {n} sites from sites_list")
        else:
            # Single-site path (unchanged — backward compatible)
            lat = site_cfg.get('latitude')
            lon = site_cfg.get('longitude')
            vs30 = site_cfg.get('vs30')
            if lat is not None and lon is not None:
                pt = Point(lon, lat)
                kwargs = {'location': pt}
                if vs30 is not None:
                    kwargs['vs30'] = vs30
                self.sitecol = SiteCollection([Site(**kwargs)])
                logger.debug(f"Initialized point site at ({lon}, {lat}) vs30={vs30}")
            else:
                raise ValueError(
                    "site_location must include latitude and longitude for hazard_curve"
                )

    def run(self):
        from openquake.fdha.calc.hazard import calculate_fdha_hazard
        return calculate_fdha_hazard(self, self.sitecol)


