# -*- coding: utf-8 -*-
"""
rank1p5_loader.py

Load Rank 1.5 rupture traces from XML (OpenQuake NRML format).

The traces are attached to the calculator as surface-like objects that expose:
  - get_fault_trace() -> Line
  - get_min_distance(SiteCollection) -> np.ndarray[km]

XML Format (OpenQuake NRML):
```xml
<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"
      xmlns:gml="http://www.opengis.net/gml">
  <rank1p5Ruptures>
    <trace name="Sibari_splay_A" sourceId="5">
      <gml:LineString>
        <gml:posList>16.338416 39.640532 16.351064 39.653699</gml:posList>
      </gml:LineString>
    </trace>
  </rank1p5Ruptures>
</nrml>
```
"""

from __future__ import annotations

import os
import logging
from typing import List, Dict, Optional
import xml.etree.ElementTree as ET
import numpy as np

from openquake.hazardlib.geo import Point, Line
from openquake.hazardlib.site import SiteCollection

logger = logging.getLogger(__name__)


class TraceOnlySurfaceAdapter:
    """
    Minimal "surface-like" adapter that wraps a polyline at the surface (z=0)
    and provides:
      - get_fault_trace() -> Line
      - get_min_distance(sitecol) -> np.ndarray of distances [km]
    
    This interface is required by decision_tree.py for case classification.
    """
    
    def __init__(
        self,
        coords_ll: List[List[float]],
        name: str = "",
        source_id: Optional[int] = None
    ):
        """
        Initialize trace adapter.
        
        Args:
            coords_ll: List of [lon, lat] coordinate pairs
            name: Trace identifier
            source_id: Reference to fault source (optional)
        """
        self.name = name
        self.source_id = source_id
        self._trace_points = [Point(lon, lat, 0.0) for lon, lat in coords_ll]
        self._line = Line(self._trace_points)
        self._np_trace = np.array(
            [[p.longitude, p.latitude] for p in self._trace_points],
            dtype=float
        )

    def get_fault_trace(self) -> Line:
        """Return the fault trace as an OpenQuake Line."""
        return self._line

    def get_min_distance(self, sitecol: SiteCollection) -> np.ndarray:
        """
        Return minimum horizontal distance [km] from each site to the trace.
        
        Uses a simple Euclidean metric in degrees scaled by 111.32 km/deg.
        
        Args:
            sitecol: SiteCollection with sites to compute distances for
            
        Returns:
            Array of distances in km, one per site
        """
        out = []
        for site in sitecol:
            s = np.array(
                [site.location.longitude, site.location.latitude],
                dtype=float
            )
            dmin_deg = float("inf")
            
            # Find minimum distance to any segment
            for i in range(len(self._np_trace) - 1):
                p1 = self._np_trace[i]
                p2 = self._np_trace[i + 1]
                seg = p2 - p1
                seglen2 = float(np.dot(seg, seg))
                
                if seglen2 == 0.0:
                    d = float(np.linalg.norm(s - p1))
                else:
                    t = float(np.clip(np.dot(s - p1, seg) / seglen2, 0.0, 1.0))
                    proj = p1 + t * seg
                    d = float(np.linalg.norm(s - proj))
                
                dmin_deg = min(d, dmin_deg)
            
            out.append(dmin_deg * 111.32)  # Convert degrees to km
        
        return np.asarray(out, dtype=float)


def load_rank1p5_traces(
    filepath: str,
    base_dir: Optional[str] = None
) -> Dict[str, TraceOnlySurfaceAdapter]:
    """
    Load rank 1.5 rupture traces from XML file.
    
    Args:
        filepath: Path to XML file
        base_dir: Base directory for relative paths
        
    Returns:
        Dictionary mapping trace name to TraceOnlySurfaceAdapter
    """
    # Resolve path
    if base_dir and not os.path.isabs(filepath):
        filepath = os.path.join(base_dir, filepath)
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Trace file not found: {filepath}")
    
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    surfaces = {}
    
    # Find trace elements (handle with or without namespace)
    for elem in root.iter():
        tag = elem.tag.split('}')[-1]  # Remove namespace prefix
        
        if tag == 'trace':
            name = elem.get('name', elem.get('id', f'trace_{len(surfaces)}'))
            source_id_str = elem.get('sourceId', elem.get('source_id'))
            source_id = int(source_id_str) if source_id_str else None
            
            # Find coordinates in gml:posList
            coords = []
            for child in elem.iter():
                child_tag = child.tag.split('}')[-1]
                if child_tag == 'posList' and child.text:
                    coords = _parse_poslist(child.text)
                    break
            
            if coords:
                surfaces[name] = TraceOnlySurfaceAdapter(
                    coords_ll=coords,
                    name=name,
                    source_id=source_id
                )
                logger.debug("Loaded trace '%s' with %d points", name, len(coords))
            else:
                logger.warning("Trace '%s' has no coordinates, skipping", name)

    logger.info("Loaded %d traces from: %s", len(surfaces), filepath)
    return surfaces


def _parse_poslist(text: str) -> List[List[float]]:
    """
    Parse GML posList (space/newline separated lon lat pairs).
    
    Format: "lon1 lat1 lon2 lat2 ..." or multiline
    """
    if not text:
        return []
    
    # Split by whitespace (handles spaces, newlines, tabs)
    values = text.strip().split()
    values = [float(v) for v in values if v]
    
    # Group into [lon, lat] pairs
    coords = []
    for i in range(0, len(values) - 1, 2):
        coords.append([values[i], values[i + 1]])
    
    return coords


def attach_rank1p5_surfaces(calculator, config_path: Optional[str] = None):
    """
    Attach rank 1.5 surfaces to calculator instance.
    
    Reads rank1p5_traces_file from config and loads the XML file.
    
    Sets:
      - calculator.rank1p5_surfaces: List[TraceOnlySurfaceAdapter]
      - calculator.rank1p5_surface_by_name: Dict[str, TraceOnlySurfaceAdapter]
    
    Args:
        calculator: Calculator instance to attach surfaces to
        config_path: Path to configuration file (for path resolution)
    """
    # Skip if already attached
    if hasattr(calculator, 'rank1p5_surfaces') and calculator.rank1p5_surfaces:
        return
    
    # Get config
    if not hasattr(calculator, 'config') or not isinstance(calculator.config, dict):
        calculator.rank1p5_surfaces = []
        calculator.rank1p5_surface_by_name = {}
        return
    
    # Get traces file path from config
    calc_cfg = calculator.config.get('calculation', calculator.config.get('parameters', {}))
    traces_file = calc_cfg.get('rank1p5_traces_file')
    
    if not traces_file:
        logger.debug("No rank1p5_traces_file in config")
        calculator.rank1p5_surfaces = []
        calculator.rank1p5_surface_by_name = {}
        return
    
    # Determine base directory
    if config_path:
        base_dir = os.path.dirname(os.path.abspath(config_path))
    elif hasattr(calculator, 'config_path'):
        base_dir = os.path.dirname(os.path.abspath(calculator.config_path))
    else:
        base_dir = os.getcwd()
    
    # Load traces
    by_name = load_rank1p5_traces(traces_file, base_dir=base_dir)
    
    # Attach to calculator
    calculator.rank1p5_surface_by_name = by_name
    calculator.rank1p5_surfaces = list(by_name.values())
    
    logger.info("Attached %d rank 1.5 surfaces to calculator", len(by_name))