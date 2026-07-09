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
        logger.debug(f"Models: primary_surf_rup={self.primary_surf_rup_model}, primary_surf_displ={self.primary_surf_displ_model}, "
                     f"secondary_surf_rup={self.secondary_surf_rup_model}, secondary_surf_displ={self.secondary_surf_displ_model}")

    @staticmethod
    def _instantiate_model(model_cfg):
        if not model_cfg or 'type' not in model_cfg:
            return None
        model_type = model_cfg['type']
        params = model_cfg.get('parameters', {})
        
        try:
            # Look up class by name in current globals or imported modules
            model_class = globals()[model_type]
            
            # Check which parameters the model's __init__ accepts
            sig = inspect.signature(model_class.__init__)
            init_params = set(sig.parameters.keys()) - {'self'}
            
            # Only pass parameters that the constructor accepts
            constructor_params = {k: v for k, v in params.items() if k in init_params}
            
            return model_class(**constructor_params)
        except Exception as e:
            logger.warning(f"Model class '{model_type}' not found or failed to instantiate: {e}")
            return None

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
    
    def _initialize_calculation_params(self):
        """Initialize calculation parameters needed by calculate_fdha_hazard."""
        # Target displacements
        target_disp = self.config.get('parameters', {}).get('target_displacement', [0.001, 0.01, 0.1, 1.0, 10.0])
        self.target_displacements = np.array(target_disp, dtype=np.float64)
        
        # Reduction configs
        self.p_sr_red_cfg = self.config.get('parameters', {}).get('primary_sr_reduction', {'method': 'median', 'q': 50})
        self.s_sr_red_cfg = self.config.get('parameters', {}).get('secondary_sr_reduction', self.p_sr_red_cfg)
        
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


