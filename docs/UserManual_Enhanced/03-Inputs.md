# Input Files

A PFDHA calculation requires a public **INI job file** plus XML logic-tree inputs. The INI file points to a **source-model logic tree** and an **FDHA model logic tree**; the source-model logic tree then points to one or more NRML source-model XML files.

This section provides a high-level overview of these inputs.

## 1. Source Model (NRML/XML)

The seismic source model defines the faults used in the analysis. It encodes where earthquakes can occur, how large they can be, and how often they happen.
This project uses NRML (Natural Hazard Risk Markup Language)—the XML schema adopted by the GEM/OpenQuake —so your source files remain interoperable with the OpenQuake Engine and related tools.

### 1.1 Source Typologies (OQ-aligned; supported in this project)

In OpenQuake, a source model may include several source types. In this PFDHA toolkit we currently support only:

-   **`SimpleFaultSource`**: A planar fault surface derived from a surface trace (polyline) plus dip and seismogenic depths; suited to shallow crustal faults.
-   **`CharacteristicFaultSource`**: Ruptures span (essentially) the entire mapped fault surface, following a characteristic magnitude/area representation. The enclosed surface may be a `simpleFaultGeometry` **or a `complexFaultGeometry`** (fault top/bottom edges); in both cases the exact NRML trace (the top edge) is retained for the FDHA distance metrics. A `complexFaultGeometry` additionally requires `[erf].complex_fault_mesh_spacing` in the job INI (see [Configuration](05-Configuration.md)).
-   **`ComplexFaultSource`**: A non-planar fault surface interpolated between explicit top and bottom edges (each a 3-D polyline), with ruptures floating over the surface as for `SimpleFaultSource`. Requires `[erf].complex_fault_mesh_spacing`.

!!! warning "Unsupported Source Types"
    Other OQ source types (e.g., `AreaSource`, `PointSource`, `MultiFaultSource`) are not supported at present.

!!! note "ComplexFaultSource and the fault trace"
    For a standalone `ComplexFaultSource`, the FDHA distance metrics (r, x/L)
    derive from the top row of the **resampled** surface mesh, not the exact
    NRML top edge — keep `complex_fault_mesh_spacing` fine near sites of
    interest. Where ruptures are meant to span the whole surface, prefer
    wrapping the `complexFaultGeometry` in a `characteristicFaultSource`,
    which retains the exact top-edge trace.

### 1.2 Required Data (fault-centric)

Each supported fault source must provide:

-   **Geometry**
    -   Surface trace as a lon–lat polyline (WGS84).
    -   Dip (degrees), upper and lower seismogenic depth (km).
    -   (Strike is implied by the trace; width is derived from dip and seismogenic thickness.)
    -   `ComplexFaultSource` / `complexFaultGeometry` instead provide explicit top and bottom edges as 3-D polylines (lon lat depth), ordered so the surface dips to the right of the strike direction (Aki & Richards convention).
-   **Kinematics**
    -   Rake (degrees). If using OQ’s distributions, provide a single nodal plane with probability 1.0.
-   **Occurrence model**
    -   `SimpleFaultSource` / `ComplexFaultSource`: typically a Gutenberg–Richter MFD (truncated) or other OQ MFD element, plus a magnitude–area relation (e.g., `WC1994`) and a rupture aspect ratio.
    -   `CharacteristicFaultSource`: a Characteristic MFD (single magnitude or narrow band) consistent with ruptures spanning the full fault.
-   **Metadata**
    -   `id`, `name`, `tectonicRegionType` (e.g., `Active Shallow Crust`).

### 1.3 How PFDHA Uses the Source Model (workflow-critical)

-   **Rupture Sampling**
    -   `SimpleFaultSource` / `ComplexFaultSource`: Magnitudes are drawn from the MFD; rupture dimensions are derived via the selected mag–area relation and aspect ratio; rupture planes are positioned along the trace/width (or the interpolated edge surface) per OQ conventions.
    -   `CharacteristicFaultSource`: Ruptures span the full mapped surface; the characteristic magnitude is used directly.
    -   Ruptures whose top edge lies deeper than 0.5 km are treated as non-surface-rupturing and are skipped: they contribute no fault-displacement hazard (applies to floating ruptures of all typologies, in both hazard-curve and hazard-map calculations).
-   **Distance Metrics for Displacement Models**
    -   The rupture surfaces/planes feed vectorized rupture–site distance calculators (e.g., to trace, to surface projection, to rupture plane).
    -   These distances, along with rake/dip, are routed through the PFDHA decision tree to select the appropriate principal and distributed displacement modelling and geometry-dependent combinations.
-   **Probability Aggregation**
    -   For each site and IM/displacement threshold, model-specific exceedance probabilities are combined (e.g., via decision-tree weights and `combine_probabilities`) to obtain the annual frequency of exceedance (AFE).
    -   !!! note "Logic Tree Support"
          Logic-tree execution is the current public path. The INI job points to `source_model_logic_tree_file` and singular `fdha_logic_tree_file`; the driver enumerates end branches, validates weights and model classes, and aggregates mean and fractile hazard results.
-   **Hazard Curve vs. Hazard Map**
    -   The hazard curve calculator operates on selected sites and IM/displacement levels.
    -   The hazard map calculator applies the same source/rupture logic across a grid of sites, honoring any distance filters to keep computations efficient.
    -   Both calculators consume the same NRML and must see consistent source `id`/parameters.

### 1.4 Example Snippet

```xml
<!-- A simple fault source defined in NRML -->
<nrml xmlns="http://openquake.org/xmlns/nrml/0.5">
  <sourceModel name="Example Source Model">
    <simpleFaultSource id="1" name="My Fault" tectonicRegion="Active Shallow Crust">
      <simpleFaultGeometry>
        <!-- Coordinates defining the fault trace -->
        <gml:LineString>
            <gml:posList>    
                -121.76 37.24
                -121.63 37.16
            </gml:posList>
        </gml:LineString>
        <dip>80.0</dip>
        <upperSeismoDepth>0.0</upperSeismoDepth>
        <lowerSeismoDepth>12.0</lowerSeismoDepth>
      </simpleFaultGeometry>
      <magScaleRel>WC1994</magScaleRel>
      <ruptAspectRatio>1.5</ruptAspectRatio>
      <incrementalMFD minMag="5.5" binWidth="0.1">
        <occurRates>0.01 0.008 0.006 0.004 0.002</occurRates>
      </incrementalMFD>
    </simpleFaultSource>
  </sourceModel>
</nrml>
```

### 1.5 Best-Practice Notes for PFDHA

-   **Stay OQ-conformant**: Use the exact MFD/geometry tags expected by the OpenQuake Engine.
-   **Keep Units Explicit**: lon/lat in degrees, depths in km, dip/rake in degrees.
-   **One Source, One Role**: Avoid mixing unsupported source types; keep your NRML to `SimpleFaultSource`, `CharacteristicFaultSource` (with simple or complex geometry), and/or `ComplexFaultSource`.
-   **Geometry Quality**: Use dense, order-consistent traces; check for self-intersections and unrealistic dips/widths.
-   **Kinematics Consistency**: Ensure rake aligns with the displacement models’ assumptions (e.g., reverse/normal/strike-slip branches).


---
## 2. Configuration File (INI)

The configuration file controls every aspect of the PFDHA calculation. It tells the `fdha` tool what kind of analysis to run, where to run it, which scientific models to use, and what parameters to apply.

-   **Format**: The configuration is written in **INI format**, following OpenQuake Engine conventions for compatibility and consistency.
-   **Content**: The INI file specifies:
    -   The calculation type is automatically detected from the `[geometry]` section (presence of `region` indicates hazard map, otherwise hazard curve).
    -   The single point site of interest for a hazard curve, or the geographical region for a hazard map.
    -   The source-model logic-tree XML file.
    -   The FDHA logic-tree XML file containing the four scientific model categories and model-specific parameters.
    -   Calculation parameters (e.g., `displacement_measure_levels`).

!!! warning "Current public INI contract"
    Public INI jobs must use `[calculation].source_model_logic_tree_file` and `[calculation].fdha_logic_tree_file`. The current loader rejects legacy public `source_model_file`, `fdha_logic_tree_files`, and `[models]` / `[models.*]` sections.

### Example Snippet

```ini
[general]
description = my_calculation

[geometry]
sites = 13.1 42.8

[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1

[calculation]
displacement_measure_levels = {"FD": [0.01, 0.1, 1.0]}
fdha_logic_tree_file = hazard_curve_minimal_fdha_logic_tree.xml
source_model_logic_tree_file = hazard_curve_minimal_source_model_logic_tree.xml
```

!!! tip "Next Step"
    The structure and all available options for the configuration file are described in detail in the **[Configuration](05-Configuration.md)** chapter.

---

## 3. Source-Model Logic Tree XML

The source-model logic tree is the public way to reference source-model XML files. A minimal wrapper uses an OpenQuake-style `sourceModel` branch:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt_sm_example">
    <logicTreeBranchSet uncertaintyType="sourceModel" branchSetID="bs_sm">
      <logicTreeBranch branchID="b_sm">
        <uncertaintyModel>source_model.xml</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
```

The source-model logic tree is parsed through OpenQuake hazardlib's `SourceModelLogicTree`, so OpenQuake source-model branch weight rules, path resolution, `applyToBranches` / `applyToSources` filters, and supported source-model uncertainty types apply.

---

## 4. FDHA Model Logic Tree XML

The FDHA model logic tree selects the scientific models. The supported FDHA uncertainty types are:

- `fdhaPrimarySRModel`
- `fdhaPrimaryFDModel`
- `fdhaSecondarySRModel`
- `fdhaSecondaryFDModel`
- `fdhaCalcRSigma` — a calculation-parameter uncertainty: each branch's `<uncertaintyModel>` carries an alternative value of the `r_sigma_km` mapping-accuracy sigma (in km, `0` allowed) rather than a model class. A job must choose one mechanism: either the scalar `[calculation].r_sigma_km` in the INI or `fdhaCalcRSigma` branch set(s) — defining both is a configuration error. See [Configuration](05-Configuration.md) for details.

!!! note "primary = principal, secondary = distributed"
    The `Primary*` uncertainty types model **principal** rupture and
    displacement on the main seismogenic fault; the `Secondary*` types model
    **distributed** (off-fault) tectonic rupture and displacement. This follows
    the IAEA primary/secondary convention; see the terminology note in
    [Scientific Models](06-Models.md) for the relationship to the
    principal/distributed terminology used by Valentini et al. (2025).

Each `<uncertaintyModel>` is an INI-like block, written as plain element text
in the same style as the OpenQuake engine's GMPE logic trees. The bracketed
line is the model class name; following lines are parameters passed to that
model. Line indentation is ignored, and wrapping the block in `<![CDATA[...]]>`
is also accepted (it parses identically).

```xml
<logicTreeBranchSet branchSetID="bs_1_primary_surf_rup"
                    uncertaintyType="fdhaPrimarySRModel">
  <logicTreeBranch branchID="B1_PRIMARY_SURF_RUP">
    <uncertaintyModel>
      [Youngs2003PrimarySR]
      style = all
    </uncertaintyModel>
    <uncertaintyWeight>1.0</uncertaintyWeight>
  </logicTreeBranch>
</logicTreeBranchSet>
```

Weights must sum to 1.0 within each branch set. The validator also checks registered model classes, valid `applyToBranches`, valid `applyToSources` when source IDs are known, and `applyToStyle` values limited to `strike-slip`, `reverse`, or `normal`.
