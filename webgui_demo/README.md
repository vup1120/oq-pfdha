# oq-pfdha — web GUI (prototype front-end)

A working web interface for
[oq-pfdha](https://github.com/vup1120/oq-pfdha), the open-source
Probabilistic Fault Displacement Hazard Analysis framework. The GUI
assembles a complete job (INI + NRML source-model logic tree + NRML FDHA
logic tree) from form inputs and executes the **actual engine**
(`FdhaLogicTree.from_ini(...).run(...)`) — the same code path as the
`fdha` command line.

> **Prototype interface** — review configurations and verify results
> independently before use in production hazard assessment.

## Features

1. **Configure**
   - Built-in source models (Norcia Case 3 / minimal example) **or upload
     your own NRML source model** — uploaded files are parsed and their
     fault sources listed (id, name, type, rake → style).
   - Hazard curve (sites) or hazard map (region + grid) geometry, fully
     editable.
   - All job parameters with the exact INI names and codebase defaults
     (cited in `app.py` comments): `investigation_time`, `r_threshold_km`,
     `near_far_threshold_km`, `rupture_mesh_spacing`, `width_of_mfd_bin`,
     `reference_vs30_value`, `return_period` (maps),
     `displacement_measure_levels` (JSON, validated live).
   - **Logic-tree builder** over the real model registry (`registry.json`):
     per-model branch weights **and editable model parameters**
     (`key = value`, the engine's native `<uncertaintyModel>` syntax),
     with the documented parameter table from the User Manual displayed
     next to each model. Live validation by the engine's own
     `nrml_reader` + `validate_spec` (FDLT rules), plus per-branch-set
     weight checks. Engine-coupled model families (e.g. the Visini et al.
     2025 distributed SR ⇔ FD models, which must co-occur in a branch) are
     chained automatically with `applyToBranches` during serialization —
     the four slots stay independent in the UI — and an **end-branch
     preview** (computed with the engine's own enumerator) shows exactly
     which model combinations and weights will run.
2. **Run** — executes the engine on a self-contained job directory under
   `webgui_demo/runs/run_<timestamp>/`; shows the generated `job.ini`,
   elapsed time, and the engine log. Errors are reported with the full
   traceback.
3. **Results** — reads the run's actual outputs: hazard curves (weighted
   mean + fractile band + optional per-branch spaghetti), displacement
   hazard maps (mean + fractile layers), validator reports, manifest, CSV
   downloads, and a ZIP of the complete output directory.

The **default configuration is a light single-branch calculation**
(Youngs 2003 chain on the Norcia Case 3 source, < 1 s for curves,
~15 s for a coarse map) — ideal for a live demonstration. Add models /
branches / sites to scale up to full studies.

## Run locally

Use a fresh virtual environment (the engine pins its own numpy/pandas
stack — do not mix with older installs):

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .                              # the oq-pfdha engine
pip install -r webgui_demo/requirements.txt   # GUI libraries
streamlit run webgui_demo/app.py
```

### Regenerating the registry snapshot

`registry.json` holds the model list, documented parameter tables and
prefills shown in the GUI. After changing the model library or the User
Manual, refresh it:

```bash
python webgui_demo/generate_registry.py
```

## Deploy on Hugging Face Spaces (Docker)

The Dockerfile must be built from the **repository root** (the image
installs the full engine):

1. Create a Space → SDK: **Docker**.
2. Push the whole repository to the Space, with `webgui_demo/Dockerfile`
   copied to the repo root (or set the Space's `dockerfile` path).
3. Spaces serves the app on port 7860 automatically.

Note: the image includes the OpenQuake engine stack (~2 GB). Runs execute
inside the Space container — size CPU accordingly and keep demo
configurations light.

## License & citation

- License: **GNU AGPL v3.0 or later** (same as oq-pfdha).
- Please cite: Chen, Y.-S. (2025). *openquake.fdha: Python tools for
  probabilistic fault displacement hazard analysis* (v1.0.0) [Software].
  Istituto Nazionale di Oceanografia e di Geofisica Sperimentale (OGS).
  https://github.com/vup1120/oq-pfdha (see `CITATION.cff`).
