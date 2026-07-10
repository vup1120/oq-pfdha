"""
Module for parsing OpenQuake NRML source models.
"""

import os
import tempfile

import numpy as np
from typing import List, Dict, Any, Optional
from openquake.commonlib.readinput import read_source_models
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import SourceModel


def parse_source_model_faults(
    fnames: Any,
    hdf5path: str = '',
    **converterparams
) -> Dict[str, SourceModel]:
    """
    Parse one or more NRML source model files into a dict of fault source objects.

    Internally calls OQ‑Engine's read_source_models, which handles
    multifault HDF5, geometry fixing, etc.

    multiFaultSource conveniences: when the caller does not pass an
    ``investigation_time``, it is sniffed from the ``<sourceModel>`` /
    ``<geometryModel>`` NRML headers (the engine SourceConverter default of
    50 yr would collide with the declared value and be rejected by the nrml
    consistency check); when ``hdf5path`` is empty but a multi-fault /
    geometry-model input is present, a temporary auxiliary .hdf5 is created
    automatically (the engine needs it to store the sections).

    :param fnames: Single filename or list of filenames of NRML XML source models.
    :param hdf5path: Path to auxiliary .hdf5 file for multifault sources.
    :param converterparams: Params passed to the converter (e.g. rupture_mesh_spacing).
    :returns: Dict mapping source_id → SourceModel instance.
    """
    # Ensure we have a list of filenames)
    if isinstance(fnames, str):
        files = [fnames]
    else:
        files = list(fnames)

    if 'investigation_time' not in converterparams:
        sniffed = _sniff_investigation_time(files)
        if sniffed is not None:
            converterparams['investigation_time'] = sniffed
    if not hdf5path and _has_multifault(files):
        # Fresh non-existing path inside a temp dir: the engine creates the
        # file itself (an empty pre-created file would fail h5py signature
        # checks).
        hdf5path = os.path.join(
            tempfile.mkdtemp(prefix='fdha_sections_'), 'sections.hdf5')

    # Read models via OQ‑Engine
    smodels = read_source_models(files, hdf5path=hdf5path, **converterparams)

    # Flatten all src_groups into a dict keyed by source_id
    sources: Dict[str, SourceModel] = {}
    for smodel in smodels:
        for sg in smodel.src_groups:
            for src in sg:
                # source_id attribute on SourceModel
                sid = getattr(src, 'source_id', None) or getattr(src, 'id', None)
                if sid is None:
                    continue
                sources[sid] = src

    _attach_original_traces(files, sources)
    return sources


def _attach_original_traces(files: List[str], sources: Dict[str, Any]) -> None:
    """Attach the exact NRML fault trace to prebuilt fault surfaces.

    The FDHA site metrics (r, rx sign, x/L) are trace quantities, but the
    hazardlib surface mesh resamples the trace at ``rupture_mesh_spacing``,
    cutting corners by up to hundreds of meters on wiggly traces (hazardlib
    even runs ``keep_corners(1.0)`` on its ``tor`` line). That error is
    irrelevant for GMPE distances but fatal for near-fault displacement
    probabilities. Retaining the raw trace here makes the FDHA distances
    exact and independent of the mesh spacing — the same decoupling the
    engine itself uses when it computes rx/ry0 from ``tor`` lines instead
    of the mesh.

    Only sources whose ruptures always span the full geometry get the trace
    (characteristicFaultSource with a prebuilt single surface). Floating
    ruptures (simpleFaultSource sub-ruptures) keep the mesh-derived trace,
    which correctly tracks the rupture extent. multiFaultSource is handled
    separately via the ECS/tor path.
    """
    for fname in files:
        try:
            root = nrml.read(fname)
        except Exception:
            continue
        for node in _walk(root):
            tag = node.tag.rsplit('}', 1)[-1]
            if tag != 'characteristicFaultSource':
                continue
            sid = node.attrib.get('id')
            src = sources.get(sid)
            if src is None or getattr(src, 'surface', None) is None:
                continue
            surface = src.surface
            # Only single-plane surfaces: a MultiSurface has no single trace.
            if hasattr(surface, 'surfaces'):
                continue
            trace = _extract_trace_coords(node)
            if trace is not None and trace.shape[0] >= 2:
                surface.original_trace = trace


def _walk(node):
    yield node
    for child in node:
        yield from _walk(child)


def _sniff_investigation_time(files: List[str]) -> Optional[float]:
    """Return the investigation_time declared in the NRML headers, if any.

    multiFaultSource / geometryModel files carry the PMF time span as an
    attribute of their root <sourceModel>/<geometryModel> node; the engine
    validates it against the converter setting, so it must be forwarded
    when the caller did not choose one explicitly.
    """
    for fname in files:
        try:
            root = nrml.read(fname)
        except Exception:
            continue
        for node in _walk(root):
            tag = node.tag.rsplit('}', 1)[-1]
            if tag in ('sourceModel', 'geometryModel'):
                value = node.attrib.get('investigation_time')
                if value is not None:
                    try:
                        return float(value)
                    except (TypeError, ValueError):
                        continue
    return None


def _has_multifault(files: List[str]) -> bool:
    """True when any input declares a multiFaultSource or geometryModel
    (both need the auxiliary sections .hdf5)."""
    for fname in files:
        try:
            root = nrml.read(fname)
        except Exception:
            continue
        for node in _walk(root):
            tag = node.tag.rsplit('}', 1)[-1]
            if tag in ('multiFaultSource', 'geometryModel'):
                return True
    return False


def _extract_trace_coords(src_node) -> 'np.ndarray | None':
    """Return the (N, 2) lon/lat trace of a fault source NRML node.

    Reads the ``simpleFaultGeometry`` LineString posList (2 values per
    vertex), or the ``complexFaultGeometry`` faultTopEdge posList (3 values
    per vertex, depth dropped).
    """
    pos_node, n_cols = None, 2
    for node in _walk(src_node):
        tag = node.tag.rsplit('}', 1)[-1]
        if tag == 'simpleFaultGeometry':
            for sub in _walk(node):
                if sub.tag.rsplit('}', 1)[-1] == 'posList':
                    pos_node, n_cols = sub, 2
                    break
            break
        if tag == 'faultTopEdge':  # complexFaultGeometry top edge
            for sub in _walk(node):
                if sub.tag.rsplit('}', 1)[-1] == 'posList':
                    pos_node, n_cols = sub, 3
                    break
            break
    if pos_node is None or pos_node.text is None:
        return None
    try:
        raw = pos_node.text
        # hazardlib's nrml reader may deliver posList pre-parsed as floats
        tokens = raw.split() if isinstance(raw, str) else raw
        vals = np.array([float(v) for v in tokens], dtype=float)
    except (ValueError, TypeError):
        return None
    if vals.size < 2 * n_cols or vals.size % n_cols != 0:
        return None
    return vals.reshape(-1, n_cols)[:, :2]
