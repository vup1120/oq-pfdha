"""
Source-model logic-tree support for FDHA.

Parsing and realisation enumeration of the logic-tree XML are delegated to
``openquake.hazardlib.logictree.SourceModelLogicTree`` so the FDHA framework
honours the same NRML schema, weight rules, multi-file branch syntax, and
modification semantics as the OpenQuake engine.  Every uncertainty type that
OpenQuake recognises in a source-model logic tree is therefore supported here:

* ``sourceModel`` (alternative full XML files);
* ``bGRRelative``, ``bGRAbsolute``, ``abGRAbsolute``, ``maxMagGRRelative``,
  ``maxMagGRRelativeNoMoBalance``, ``maxMagGRAbsolute``, ``abMaxMagAbsolute``,
  ``incrementalMFDAbsolute``, ``truncatedGRFromSlipAbsolute``;
* ``simpleFaultGeometryAbsolute``, ``simpleFaultDipRelative``,
  ``simpleFaultDipAbsolute``, ``complexFaultGeometryAbsolute``,
  ``characteristicFaultGeometryAbsolute``, ``areaSourceGeometryAbsolute``;
* ``setMSRAbsolute``, ``setLowerSeismDepthAbsolute``,
  ``setUpperSeismDepthAbsolute``, ``recomputeMmax``, ``dummy``.

Each realisation yielded by OpenQuake's ``SourceModelLogicTree`` iterator
carries the base source-model XML file(s) chosen by the realisation's
``sourceModel`` branch plus a list of
``(uncertainty_type, value, applyToSources)`` modifications.
The driver applies them via :func:`apply_realization_to_sources` (a thin
wrapper around :func:`openquake.hazardlib.lt.apply_uncertainties`) before
handing the resulting source dict to the FDHA calculator.

A backward-compatible :func:`load_source_model_branches` returns a single
implicit ``sm0`` realisation when only ``source_model_file`` is configured,
so callers see a uniform list of :class:`SourceModelBranch` regardless of
whether a logic tree is in use.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

# Conceptual branch ID assigned to the implicit branch built from a single
# ``source_model_file`` entry (no logic tree).
IMPLICIT_BRANCH_ID = "sm0"

# Tolerance kept for backward compatibility with callers; OpenQuake itself
# enforces weight sums internally during parsing.
WEIGHT_SUM_TOLERANCE = 1e-6


class SourceModelLogicTreeError(Exception):
    """Raised for any parsing/validation problem in a source-model logic tree.

    Wraps :class:`openquake.hazardlib.logictree.LogicTreeError` and
    :class:`openquake.baselib.InvalidFile` so callers can catch a single,
    FDHA-specific exception type without depending on hazardlib internals.
    """


@dataclass(frozen=True)
class UncertaintyApplication:
    """One NRML modification to apply to a source group.

    ``apply_to_sources`` and ``apply_to_branches`` reproduce the OpenQuake
    branch-set filters that decided whether this branch participates in the
    current realisation; they are kept here for traceability/manifests.
    """

    uncertainty_type: str
    value: Any
    branch_set_id: str
    branch_id: str
    apply_to_sources: Optional[tuple[str, ...]] = None
    apply_to_branches: Optional[tuple[str, ...]] = None


@dataclass(frozen=True)
class SourceModelBranch:
    """One realisation through the source-model logic tree.

    The name is kept for backward compatibility with the earlier
    sourceModel-only design; semantically each instance is now a *full*
    realisation that combines:

    * a base source-model XML file (or whitespace-separated set of files)
      coming from the realisation's ``sourceModel`` branch, and
    * an ordered list of NRML uncertainty modifications to apply on top.

    ``branch_id`` is composed by joining each branch ID along the realisation
    path with ``"|"`` so it remains stable across runs.  ``weight`` is the
    product of all branch weights along the path (sum across realisations is
    1.0 within OpenQuake's tolerance).
    """

    branch_id: str
    source_model_file: str
    weight: float
    metadata: dict = field(default_factory=dict)
    uncertainties: tuple[UncertaintyApplication, ...] = ()
    # Position of this realisation in the enumerated set, mirroring OpenQuake's
    # ``SourceModelLogicTree`` realisation ordinals.  Default ``0`` keeps the
    # implicit (single-branch) case OpenQuake-compatible without forcing
    # callers to set it.
    ordinal: int = 0


# ---------------------------------------------------------------- parsing


def parse_source_model_logic_tree(
    xml_path: Union[str, Path],
) -> list[SourceModelBranch]:
    """Parse a source-model logic-tree and return one realisation per path.

    Internally delegates to
    :class:`openquake.hazardlib.logictree.SourceModelLogicTree`, whose
    iterator enumerates (or samples) realisations, and to its
    :meth:`bset_values` for the per-realisation ``(branchset, value)`` pairs.
    The whole OpenQuake NRML schema therefore applies transparently: weight
    sums, relative-path/file-existence validation (raised at construction
    time), ``applyToBranches``, ``applyToSources``, and every supported
    ``uncertaintyType``.  We only translate OpenQuake's ``Realization`` into
    the FDHA-facing :class:`SourceModelBranch` record and resolve the
    sourceModel branch's XML file(s) to absolute paths for the manifests.
    """
    xml_path = Path(xml_path)
    if not xml_path.exists():
        raise SourceModelLogicTreeError(
            f"Source-model logic tree file not found: {xml_path}"
        )

    smlt = _open_smlt(xml_path)
    basedir = Path(smlt.basepath)

    realizations: list[SourceModelBranch] = []
    for rlz in smlt:
        path_ids = list(rlz.lt_path)
        sm_branch_id = path_ids[0]
        # ``rlz.value`` is one entry per branchset; the first is the
        # sourceModel branch's ``<uncertaintyModel>`` (a whitespace-separated
        # set of one or more XML files, already validated by hazardlib).
        files_raw = _split_branch_value(rlz.value[0])
        abs_files = [str((basedir / f).resolve()) for f in files_raw]
        joined = abs_files[0] if len(abs_files) == 1 else ",".join(abs_files)

        uncertainties: list[UncertaintyApplication] = []
        for bset, value in smlt.bset_values(path_ids):
            ordinal = getattr(bset, "ordinal", None)
            br_id = path_ids[ordinal] if ordinal is not None else ""
            filters = getattr(bset, "filters", {})
            uncertainties.append(
                UncertaintyApplication(
                    uncertainty_type=bset.uncertainty_type,
                    value=value,
                    branch_set_id=getattr(bset, "id", "") or "",
                    branch_id=br_id,
                    apply_to_sources=_tuple_or_none(
                        filters.get("applyToSources")
                    ),
                    apply_to_branches=_tuple_or_none(
                        filters.get("applyToBranches")
                    ),
                )
            )

        realizations.append(
            SourceModelBranch(
                branch_id="|".join(path_ids),
                source_model_file=joined,
                weight=float(rlz.weight),
                metadata={
                    "branch_path": path_ids,
                    "files": abs_files,
                    "logic_tree_file": str(xml_path.resolve()),
                    "source_model_branch_id": sm_branch_id,
                    "raw_value": rlz.value[0],
                },
                uncertainties=tuple(uncertainties),
                ordinal=rlz.ordinal,
            )
        )

    if not realizations:
        raise SourceModelLogicTreeError(
            f"No realisations enumerated from {xml_path}"
        )

    return realizations


def _open_smlt(xml_path: Path):
    try:
        from openquake.baselib import InvalidFile
        from openquake.hazardlib.logictree import (
            LogicTreeError,
            SourceModelLogicTree,
        )
    except ImportError as e:  # pragma: no cover
        raise SourceModelLogicTreeError(
            f"openquake.hazardlib is required to parse source-model logic "
            f"trees but failed to import: {e}"
        )
    try:
        return SourceModelLogicTree(str(xml_path))
    except (LogicTreeError, InvalidFile) as e:
        raise SourceModelLogicTreeError(str(e)) from e


def _split_branch_value(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value if str(v).strip()]
    text = str(value).strip()
    return text.split() if text else []


def _tuple_or_none(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple, set)):
        return tuple(str(v) for v in value)
    return (str(value),)


# -------------------------------------------------------- application API


def apply_realization_to_sources(
    branch: SourceModelBranch,
    fault_sources: dict[str, Any],
    *,
    converter_params: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Return a copy of ``fault_sources`` with this branch's NRML
    modifications applied.

    For every :class:`UncertaintyApplication` recorded on the realisation we
    build a one-element OpenQuake ``BranchSet`` so we can reuse the engine's
    canonical :func:`apply_uncertainties` machinery (which honours
    ``applyToSources`` and source-type filters per uncertainty type).

    Sources that are not modified are returned unchanged (by reference);
    sources that are modified are deep-copied before the modification so the
    original calculator state is never mutated.
    """
    if not branch.uncertainties:
        return dict(fault_sources)

    try:
        import copy as _copy

        from openquake.hazardlib import lt as _lt
    except ImportError as e:  # pragma: no cover
        raise SourceModelLogicTreeError(
            f"openquake.hazardlib is required to apply uncertainties: {e}"
        )

    out: dict[str, Any] = dict(fault_sources)
    for ua in branch.uncertainties:
        try:
            bset = _build_pseudo_branchset(ua)
        except NotImplementedError as e:
            # hazardlib's BranchSet rejects unrecognised uncertainty types at
            # construction time; keep the documented FDHA exception type.
            raise SourceModelLogicTreeError(
                f"Unknown uncertainty type '{ua.uncertainty_type}' "
                f"in branch '{ua.branch_id}'"
            ) from e
        for sid, src in list(out.items()):
            if not bset.filter_source(src):
                continue
            new_src = _copy.deepcopy(src)
            try:
                _lt.apply_uncertainty(ua.uncertainty_type, new_src, ua.value)
            except KeyError:
                raise SourceModelLogicTreeError(
                    f"Unknown uncertainty type '{ua.uncertainty_type}' "
                    f"in branch '{ua.branch_id}'"
                )
            out[sid] = new_src
    return out


def _build_pseudo_branchset(ua: UncertaintyApplication):
    """Build a minimal OpenQuake ``BranchSet`` whose ``filter_source`` honours
    the realisation's ``applyToSources`` filter.

    We only need ``filter_source`` to behave correctly; we do not actually
    iterate this BranchSet's branches (the value to apply is already known
    from ``UncertaintyApplication``).
    """
    from openquake.hazardlib import lt as _lt

    filters: dict[str, Any] = {}
    if ua.apply_to_sources:
        filters["applyToSources"] = list(ua.apply_to_sources)
    if ua.apply_to_branches:
        filters["applyToBranches"] = list(ua.apply_to_branches)
    bset = _lt.BranchSet(ua.uncertainty_type, filters=filters)
    return bset


# ----------------------------------------------------------- normalisation


def load_source_model_branches(
    config: dict[str, Any],
    config_dir: Union[str, Path],
) -> list[SourceModelBranch]:
    """Return the source-model realisations implied by an FDHA configuration.

    - When ``[calculation].source_model_logic_tree_file`` is set, parse it via
      OpenQuake's :class:`SourceModelLogicTree` and return one realisation per
      enumerated path.
    - Otherwise fall back to ``[calculation].source_model_file`` (or
      ``[calculation].source_model``) and synthesise a single implicit
      realisation with id ``"sm0"`` and weight ``1.0`` so downstream code can
      treat the two cases uniformly.
    """
    config_dir = Path(config_dir)
    calc_cfg = config.get("calculation", config.get("parameters", {})) or {}

    smlt_file = calc_cfg.get("source_model_logic_tree_file")
    if smlt_file:
        smlt_path = Path(smlt_file)
        if not smlt_path.is_absolute():
            smlt_path = config_dir / smlt_path
        return parse_source_model_logic_tree(smlt_path)

    sm = calc_cfg.get("source_model_file", calc_cfg.get("source_model"))
    if not sm:
        raise SourceModelLogicTreeError(
            "No source model specified: set 'source_model_file' or "
            "'source_model_logic_tree_file' under [calculation]."
        )

    if isinstance(sm, (list, tuple)):
        sm_paths = [_resolve_or_raise(entry, config_dir) for entry in sm]
        joined = ",".join(sm_paths)
        return [
            SourceModelBranch(
                branch_id=IMPLICIT_BRANCH_ID,
                source_model_file=joined,
                weight=1.0,
                metadata={"implicit": True, "files": sm_paths},
            )
        ]

    resolved = _resolve_or_raise(sm, config_dir)
    return [
        SourceModelBranch(
            branch_id=IMPLICIT_BRANCH_ID,
            source_model_file=resolved,
            weight=1.0,
            metadata={"implicit": True, "files": [resolved]},
        )
    ]


def _resolve_or_raise(raw: str, config_dir: Path) -> str:
    p = Path(raw)
    if not p.is_absolute():
        candidate = config_dir / p
        if candidate.exists():
            p = candidate
    if not p.exists():
        raise SourceModelLogicTreeError(
            f"Source model file not found: {raw} (looked in {config_dir})"
        )
    return str(p.resolve())


def expand_branch_paths(branch: SourceModelBranch) -> list[str]:
    """Return the source-model XML files implied by a branch."""
    files = branch.metadata.get("files") if branch.metadata else None
    if files:
        return list(files)
    if "," in branch.source_model_file:
        return [p for p in branch.source_model_file.split(",") if p]
    return [branch.source_model_file]


# ----------------------------------------------------- OpenQuake-style API


def prepare_source_model_realizations(
    config: dict[str, Any],
    config_dir: Union[str, Path],
) -> list[SourceModelBranch]:
    """OpenQuake-style entry point: return one source-model realisation per
    enumerated path, regardless of whether a logic tree is configured.

    Mirrors :func:`openquake.hazardlib.get_smlt`/
    :meth:`SourceModelLogicTree.trivial`: when only ``source_model_file`` is
    set we return a single trivial realisation (ordinal 0, weight 1.0); when
    ``source_model_logic_tree_file`` is set we delegate parsing to OpenQuake
    and return one ``SourceModelBranch`` per realisation, each with its
    branch path, ordinal, and weight.

    The returned objects are accepted by :func:`apply_realization_to_sources`
    and by both the hazard-curve and hazard-map driver paths, so callers can
    write a single ``for sm_real in prepare_source_model_realizations(...)``
    loop and have map mode honour SMLT exactly the way curve mode does.
    """
    branches = load_source_model_branches(config, config_dir)
    return [
        SourceModelBranch(
            branch_id=b.branch_id,
            source_model_file=b.source_model_file,
            weight=b.weight,
            metadata=b.metadata,
            uncertainties=b.uncertainties,
            ordinal=i,
        )
        for i, b in enumerate(branches)
    ]
