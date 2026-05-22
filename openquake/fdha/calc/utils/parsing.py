"""
Module for parsing OpenQuake NRML source models.
"""

from typing import List, Dict, Any
from openquake.commonlib.readinput import read_source_models
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
    return sources
