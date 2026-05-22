"""
Configuration loader for FDHA calculations

This module provides unified configuration loading and validation functionality
for fault displacement hazard analysis calculations.
The public v5 configuration surface is INI-only.
"""

import json
import configparser
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
from dataclasses import dataclass, field


# Canonical example for OpenQuake-style FDHA classical + logic-tree layout.
_FDHA_CANONICAL_JOB_INI_REL = Path(
    "examples/hazard_curve_minimal.ini"
)


class ConfigurationError(Exception):
    """Base exception for configuration errors"""
    pass


class ConfigValidationError(ConfigurationError):
    """Exception raised when configuration validation fails"""
    pass


@dataclass
class ERFConfig:
    """Earthquake Rupture Forecast configuration"""
    rupture_mesh_spacing: float = 0.5
    width_of_mfd_bin: float = 0.1
    
    def __post_init__(self):
        """Validate configuration after initialization"""
        if self.rupture_mesh_spacing <= 0:
            raise ConfigValidationError("rupture_mesh_spacing must be > 0")
        if self.width_of_mfd_bin <= 0:
            raise ConfigValidationError("width_of_mfd_bin must be > 0")


@dataclass
class CalculationConfig:
    """Calculation configuration"""
    target_displacements: list = field(default_factory=lambda: [0.001, 0.01, 0.1, 1.0, 10.0])
    investigation_time: float = 1.0
    rupture_mesh_spacing: float = 2.0
    source_model_file: Optional[str] = None
    source_model_logic_tree_file: Optional[str] = None

    def get_source_model_paths(self) -> List[str]:
        """
        Get list of source model paths from configuration.
        Returns a list containing source_model_file or source_model_logic_tree_file if present.
        """
        paths = []
        if self.source_model_file:
            paths.append(self.source_model_file)
        if self.source_model_logic_tree_file:
            paths.append(self.source_model_logic_tree_file)
        return paths
    
    def __post_init__(self):
        """Validate configuration after initialization"""
        if self.investigation_time <= 0:
            raise ConfigValidationError("investigation_time must be > 0")
        if self.rupture_mesh_spacing <= 0:
            raise ConfigValidationError("rupture_mesh_spacing must be > 0")
        if not self.target_displacements:
            raise ConfigValidationError("target_displacements cannot be empty")


@dataclass
class FDHAConfiguration:
    """Main FDHA configuration container"""
    erf: ERFConfig = field(default_factory=ERFConfig)
    calculation: CalculationConfig = field(default_factory=CalculationConfig)
    source_model_file: Optional[str] = None
    output_dir: Optional[str] = None
    
    def __post_init__(self):
        """Validate configuration after initialization"""
        if self.source_model_file and not Path(self.source_model_file).exists():
            raise ConfigValidationError(f"Source model file not found: {self.source_model_file}")


def load_config(config_file: Union[str, Path]) -> Dict[str, Any]:
    """
    Unified configuration loader for v5 canonical INI files.
    
    This is the canonical configuration loader that should be used throughout
    the codebase.
    
    Args:
        config_file: Path to configuration file (INI)
        
    Returns:
        Configuration dictionary with nested structure
        
    Raises:
        ConfigurationError: If file cannot be read or parsed
    """
    config_path = Path(config_file)
    
    if not config_path.exists():
        raise ConfigurationError(f"Configuration file not found: {config_file}")
    
    ext = config_path.suffix.lower()
    
    if ext == '.ini':
        return _load_ini_config(config_path)
    if ext == '.toml':
        raise ConfigurationError(
            "TOML configuration files are no longer supported. "
            "Use a v5 canonical INI job with `fdha_logic_tree_file` and "
            "`source_model_logic_tree_file`."
        )
    raise ConfigurationError(
        f"Unsupported configuration file extension {ext!r}; use a v5 canonical .ini job."
    )


def _load_ini_config(config_path: Path) -> Dict[str, Any]:
    """
    Load configuration from INI file with support for nested sections and JSON values.
    
    Handles:
    - Nested sections: [models.primary_surf_rup] -> config['models']['primary_surf_rup']
    - JSON values: displacement_measure_levels = {"FD": [...]}
    - Numbers, booleans, None
    - Normalization to match TOML expected format
    """
    cp = configparser.RawConfigParser()
    cp.optionxform = str  # Preserve key case
    cp.read(config_path)

    _reject_legacy_public_ini(cp, config_path)
    
    config = {}
    
    for section in cp.sections():
        # Split section name by dots for nested structure
        section_parts = section.split('.')
        
        # Navigate/create nested structure
        current = config
        for part in section_parts:
            if part not in current:
                current[part] = {}
            current = current[part]
        
        # Parse each key-value pair
        for key, value in cp.items(section):
            current[key] = _parse_ini_value(value)
    
    # Normalize INI-specific format to match TOML expected structure
    _normalize_ini_config(config, config_path)
    
    return config


def _is_materialized_branch_ini(config_path: Path) -> bool:
    """Internal logic-tree branch INIs still feed the calculator directly."""
    return "branch_configs" in config_path.parts


def _reject_legacy_public_ini(cp: configparser.RawConfigParser, config_path: Path) -> None:
    """Reject legacy user-facing INI constructs before normalisation."""
    if _is_materialized_branch_ini(config_path):
        return

    if cp.has_section("calculation"):
        if cp.has_option("calculation", "source_model_file"):
            raise ConfigurationError(
                "Legacy key `source_model_file` is no longer supported. "
                "Replace with `source_model_logic_tree_file = <wrapper.xml>`. "
                "See CHANGELOG section vX.Y."
            )
        if cp.has_option("calculation", "fdha_logic_tree_files"):
            raise ConfigurationError(
                "Legacy key `fdha_logic_tree_files` (plural) is no longer supported. "
                "Replace with `fdha_logic_tree_file = <merged.xml>` containing all "
                "styles via `applyToStyle`. See CHANGELOG section vX.Y."
            )

    for section in cp.sections():
        if section == "models" or section.startswith("models."):
            raise ConfigurationError(
                "Legacy `[models]` section is no longer supported. Model parameters "
                "are now defined inside `<uncertaintyModel>` blocks of the FDHA logic "
                "tree XML. See CHANGELOG section vX.Y."
            )

    if cp.has_section("parameters") and cp.has_option("parameters", "target_displacement"):
        raise ConfigurationError(
            "Legacy key `target_displacement` is no longer supported. Replace with "
            "displacement_measure_levels = {\"FD\": [...]}. See CHANGELOG section vX.Y."
        )


def canonical_fdha_job_ini_hint() -> str:
    """Reference path cited in ConfigurationError messages (repo-relative)."""
    return str(_FDHA_CANONICAL_JOB_INI_REL)


def calculation_requests_fdha_logic_tree(calc: Dict[str, Any]) -> bool:
    """True if ``[calculation]`` declares the canonical FDHA logic tree file."""
    if not isinstance(calc, dict):
        return False
    sg = calc.get("fdha_logic_tree_file")
    if isinstance(sg, str) and sg.strip():
        return True
    return False


def is_canonical_fdha_logic_tree_ini(calc: Dict[str, Any]) -> bool:
    """Canonical: ``source_model_logic_tree_file`` + ``fdha_logic_tree_file`` (singular only)."""
    if not isinstance(calc, dict):
        return False
    return bool(
        calc.get("source_model_logic_tree_file")
        and calc.get("fdha_logic_tree_file")
        and not calc.get("fdha_logic_tree_files")
    )


def validate_public_logic_tree_ini_file(config_path: Path, config: Dict[str, Any]) -> None:
    """Raise ConfigurationError when a user's logic-tree ``job.ini`` violates v5 layout.

    Branch materialised INIs (written under ``branch_configs/`` with ``[models.*]``)
    are not passed through this function — only the top-level job the user passes
    to ``FdhaLogicTree.from_ini`` / ``fdha job.ini``.

    Legacy user-facing keys are rejected during raw INI loading.
    """
    if config_path.suffix.lower() != ".ini":
        return

    calc = config.get("calculation", {})
    if not isinstance(calc, dict) or not calculation_requests_fdha_logic_tree(calc):
        return

    rcp = configparser.RawConfigParser()
    rcp.optionxform = str
    read_ok = rcp.read(config_path)
    if not read_ok:
        return

    canonical = is_canonical_fdha_logic_tree_ini(calc)

    if canonical:
        if not str(calc.get("fdha_logic_tree_file", "")).strip():
            raise ConfigurationError(
                "Canonical logic-tree jobs require [calculation].fdha_logic_tree_file."
            )
    else:
        raise ConfigurationError(
            "Canonical logic-tree jobs require both "
            "[calculation].source_model_logic_tree_file and "
            "[calculation].fdha_logic_tree_file. "
            f"See canonical layout in {canonical_fdha_job_ini_hint()}."
        )


def _normalize_ini_config(config: Dict[str, Any], config_path: Path) -> None:
    """
    Normalize INI configuration to match TOML expected format.
    
    Handles:
    - [geometry].corner_points (string) -> [site_location].corner_points (list)
    - [site_params].reference_vs30_value -> [site_location].vs30
    - [calculation].displacement_measure_levels -> [parameters].target_displacement
    - [calculation].rank1p5_traces_file parsing
    """
    import os
    
    # 1. Convert geometry corner_points string to site_location corner_points list
    if 'geometry' in config and 'corner_points' in config['geometry']:
        if 'site_location' not in config:
            config['site_location'] = {}
        
        corner_str = config['geometry']['corner_points']
        if isinstance(corner_str, str):
            # Parse "lon1 lat1, lon2 lat2, ..." format
            corners = []
            for pt in corner_str.split(','):
                pt = pt.strip()
                if pt:
                    parts = pt.split()
                    if len(parts) >= 2:
                        corners.append([float(parts[0]), float(parts[1])])
            if corners:
                config['site_location']['corner_points'] = corners
                # Calculate centroid and set as latitude/longitude for point site
                import numpy as np
                corners_array = np.array(corners)
                centroid_lon = float(np.mean(corners_array[:, 0]))
                centroid_lat = float(np.mean(corners_array[:, 1]))
                config['site_location']['longitude'] = centroid_lon
                config['site_location']['latitude'] = centroid_lat
    
    # 2. Convert site_params.reference_vs30_value to site_location.vs30
    if 'site_params' in config and 'reference_vs30_value' in config['site_params']:
        if 'site_location' not in config:
            config['site_location'] = {}
        if 'vs30' not in config['site_location']:
            config['site_location']['vs30'] = config['site_params']['reference_vs30_value']
    
    # 3. Convert calculation.displacement_measure_levels to parameters.target_displacement
    if 'calculation' in config and 'displacement_measure_levels' in config['calculation']:
        if 'parameters' not in config:
            config['parameters'] = {}
        
        dml = config['calculation']['displacement_measure_levels']
        if isinstance(dml, dict) and 'FD' in dml:
            if 'target_displacement' not in config['parameters']:
                config['parameters']['target_displacement'] = dml['FD']

    # 3b. OpenQuake-style [calculation].max_distance_km → [geometry].max_distance_km
    if 'calculation' in config and 'max_distance_km' in config['calculation']:
        if 'geometry' not in config:
            config['geometry'] = {}
        if 'max_distance_km' not in config['geometry']:
            config['geometry']['max_distance_km'] = config['calculation']['max_distance_km']
    
    # 4. Copy case from parameters or calculation section
    if 'parameters' in config and 'case' in config['parameters']:
        pass  # Already set
    elif 'calculation' in config and 'case' in config['calculation']:
        if 'parameters' not in config:
            config['parameters'] = {}
        config['parameters']['case'] = config['calculation']['case']
    
    # 5. Handle rank1p5_traces_file if present
    if 'calculation' in config and 'rank1p5_traces_file' in config['calculation']:
        traces_file = config['calculation']['rank1p5_traces_file']
        if traces_file and 'rank1p5_ruptures' not in config:
            # Parse the XML file for rank1.5 traces
            traces_path = config_path.parent / traces_file
            if traces_path.exists():
                try:
                    traces = _parse_rank1p5_traces_xml(str(traces_path))
                    if traces:
                        config['rank1p5_ruptures'] = {'trace': traces}
                except Exception:
                    pass  # Silently ignore parsing errors
    
    # 6. Copy near_far_threshold_km from calculation to parameters
    if 'calculation' in config and 'near_far_threshold_km' in config['calculation']:
        if 'parameters' not in config:
            config['parameters'] = {}
        if 'near_far_threshold_km' not in config['parameters']:
            config['parameters']['near_far_threshold_km'] = config['calculation']['near_far_threshold_km']
    
    # 7. Convert geometry.sites string to site_location.latitude/longitude
    #    Supports both single-site "lon lat" and multi-site "lon1 lat1, lon2 lat2 [,depth]..."
    if 'geometry' in config and 'sites' in config['geometry']:
        if 'site_location' not in config:
            config['site_location'] = {}

        sites_str = config['geometry']['sites']
        if isinstance(sites_str, str):
            sites_list = _parse_sites_string(sites_str)
            if sites_list:
                config['site_location']['sites_list'] = sites_list
                # Backward compat: always expose first site as scalar lon/lat
                config['site_location']['longitude'] = sites_list[0]['longitude']
                config['site_location']['latitude'] = sites_list[0]['latitude']

    # 8. Convert geometry.sites_csv path to site_location.sites_list
    if 'geometry' in config and 'sites_csv' in config['geometry']:
        if 'site_location' not in config:
            config['site_location'] = {}

        csv_val = config['geometry']['sites_csv']
        if isinstance(csv_val, str) and csv_val.strip():
            csv_path = config_path.parent / csv_val.strip()
            sites_list = _parse_sites_csv(csv_path)
            if sites_list:
                config['site_location']['sites_list'] = sites_list
                config['site_location']['longitude'] = sites_list[0]['longitude']
                config['site_location']['latitude'] = sites_list[0]['latitude']


def _parse_rank1p5_traces_xml(xml_path: str) -> List[Dict[str, Any]]:
    """Parse rank1p5 traces from XML file."""
    import xml.etree.ElementTree as ET
    
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Handle namespaces - the elements have full namespace URLs
    nrml_ns = 'http://openquake.org/xmlns/nrml/0.5'
    gml_ns = 'http://www.opengis.net/gml'
    
    traces = []
    
    # Find all trace elements under rank1p5Ruptures
    # Use iteration with namespace-qualified tags
    for trace in root.iter(f'{{{nrml_ns}}}trace'):
        name = trace.get('name')
        # Find posList within LineString
        pos_list = trace.find(f'.//{{{gml_ns}}}posList')
        if pos_list is not None and pos_list.text:
            coords_text = pos_list.text.strip().split()
            # coords are "lon1 lat1 lon2 lat2 ..."
            coords = []
            for i in range(0, len(coords_text), 2):
                if i + 1 < len(coords_text):
                    coords.append([float(coords_text[i]), float(coords_text[i + 1])])
            if name and coords:
                traces.append({
                    'name': name,
                    'geometry': {'type': 'Line', 'coords': coords}
                })
    
    # Try without namespace if no traces found
    if not traces:
        for trace in root.iter('trace'):
            name = trace.get('name')
            pos_list = trace.find(f'.//{{{gml_ns}}}posList')
            if pos_list is None:
                pos_list = trace.find('.//posList')
            if pos_list is not None and pos_list.text:
                coords_text = pos_list.text.strip().split()
                coords = []
                for i in range(0, len(coords_text), 2):
                    if i + 1 < len(coords_text):
                        coords.append([float(coords_text[i]), float(coords_text[i + 1])])
                if name and coords:
                    traces.append({
                        'name': name,
                        'geometry': {'type': 'Line', 'coords': coords}
                    })
    
    return traces


def _validate_lonlat(lon: float, lat: float) -> None:
    """Raise ConfigValidationError if coordinates are out of range."""
    if not (-180.0 < lon <= 360.0):
        raise ConfigValidationError(
            f"Longitude {lon} is out of range (-180, 360]. Check [geometry].sites or sites_csv."
        )
    if not (-90.0 <= lat <= 90.0):
        raise ConfigValidationError(
            f"Latitude {lat} is out of range [-90, 90]. Check [geometry].sites or sites_csv."
        )


def _parse_sites_string(sites_str: str) -> List[Dict[str, Any]]:
    """
    Parse OpenQuake-style multi-site string: "lon1 lat1 [depth1], lon2 lat2 [depth2], ..."

    Each comma-separated token is one site; whitespace within a token separates
    lon, lat, and optional depth fields.  Single-site "lon lat" is handled as a
    degenerate case (no comma required).

    Returns a list of dicts with keys 'longitude', 'latitude', and optionally 'depth'.
    """
    sites_list = []
    for token in sites_str.split(','):
        token = token.strip()
        if not token:
            continue
        parts = token.split()
        if len(parts) < 2:
            raise ConfigValidationError(
                f"Each site must have at least lon and lat; got: {token!r}"
            )
        lon = float(parts[0])
        lat = float(parts[1])
        _validate_lonlat(lon, lat)
        entry: Dict[str, Any] = {'longitude': lon, 'latitude': lat}
        if len(parts) >= 3:
            entry['depth'] = float(parts[2])
        sites_list.append(entry)
    return sites_list


def _parse_sites_csv(csv_path) -> List[Dict[str, Any]]:
    """
    Parse a CSV file with header site_id,lon,lat[,depth].

    'lon' and 'lat' columns are required; 'depth' is optional.
    Raises ConfigValidationError if the file is missing, unreadable, or lacks
    required columns.
    """
    import csv

    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise ConfigValidationError(f"sites_csv file not found: {csv_path}")

    sites_list = []
    try:
        with open(csv_path, newline='') as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                raise ConfigValidationError(f"sites_csv file is empty: {csv_path}")
            # Normalise header names to lowercase for robust matching
            fieldnames_lower = [f.strip().lower() for f in reader.fieldnames]
            if 'lon' not in fieldnames_lower or 'lat' not in fieldnames_lower:
                raise ConfigValidationError(
                    f"sites_csv must have 'lon' and 'lat' columns; "
                    f"found: {reader.fieldnames} in {csv_path}"
                )
            has_depth = 'depth' in fieldnames_lower

            for i, row in enumerate(reader):
                # Use lowercase normalised keys
                row_lower = {k.strip().lower(): v.strip() for k, v in row.items()}
                try:
                    lon = float(row_lower['lon'])
                    lat = float(row_lower['lat'])
                except (ValueError, KeyError) as exc:
                    raise ConfigValidationError(
                        f"sites_csv row {i + 2} has invalid lon/lat: {row} ({exc})"
                    )
                _validate_lonlat(lon, lat)
                entry: Dict[str, Any] = {'longitude': lon, 'latitude': lat}
                if has_depth and row_lower.get('depth', '').strip():
                    try:
                        entry['depth'] = float(row_lower['depth'])
                    except ValueError as exc:
                        raise ConfigValidationError(
                            f"sites_csv row {i + 2} has invalid depth: {row} ({exc})"
                        )
                sites_list.append(entry)
    except ConfigValidationError:
        raise
    except Exception as exc:
        raise ConfigValidationError(f"Failed to read sites_csv {csv_path}: {exc}")

    if not sites_list:
        raise ConfigValidationError(f"sites_csv file has no data rows: {csv_path}")

    return sites_list


def _parse_ini_value(value: str) -> Any:
    """
    Parse INI value, handling JSON, numbers, booleans.
    
    This matches the logic from calculators.py._parse_ini_value() to ensure
    consistent behavior.
    """
    value = value.strip()
    
    if not value:
        return value
    
    # Boolean
    if value.lower() in ('true', 'yes', 'on'):
        return True
    if value.lower() in ('false', 'no', 'off'):
        return False
    
    # None
    if value.lower() in ('none', 'null'):
        return None
    
    # Try JSON first (for lists, dicts, etc.)
    try:
        return json.loads(value.replace("'", '"'))
    except (json.JSONDecodeError, ValueError):
        pass
    
    # Try numeric conversion
    try:
        if '.' in value:
            return float(value)
        else:
            return int(value)
    except ValueError:
        pass
    
    # Strip quotes if present
    if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    elif len(value) >= 2 and value.startswith("'") and value.endswith("'"):
        value = value[1:-1]
    
    # Keep as string
    return value


def load_fdha_config(config_file: Union[str, Path]) -> FDHAConfiguration:
    """Legacy TOML configuration entry point kept only as an explicit hard stop."""
    raise ConfigurationError(
        "TOML configuration files are no longer supported. "
        "Use `load_config()` with a v5 canonical INI job instead."
    )


def resolve_path(path: Union[str, Path], base_dir: Optional[Union[str, Path]] = None) -> Path:
    """
    Resolve relative paths to absolute paths
    
    Args:
        path: Path to resolve
        base_dir: Base directory for relative paths
        
    Returns:
        Resolved absolute path
    """
    path = Path(path)
    
    if path.is_absolute():
        return path
    
    if base_dir is None:
        base_dir = Path.cwd()
    else:
        base_dir = Path(base_dir)
    
    return base_dir / path
