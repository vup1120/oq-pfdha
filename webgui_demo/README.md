# oq-pfdha — web GUI demonstration prototype

A **visual proof-of-concept** of a web interface for
[oq-pfdha](https://github.com/vup1120/oq-pfdha), the open-source
Probabilistic Fault Displacement Hazard Analysis framework built on
OpenQuake Engine infrastructure.

> ⚠️ **Prototype — not for production hazard assessment.**
> No real calculation is performed. The "Run" button simulates progress and
> the Results page displays **pre-computed** Norcia Case 3 outputs
> (IAEA TECDOC-2092 benchmark) committed under `data/norcia_case3/`
> (provenance: `data/norcia_case3/PROVENANCE.md`). Changing models or
> parameters does **not** change the displayed results.

## What it demonstrates

1. **Configure** — source-model selection, calculation parameters with the
   exact INI names and codebase defaults (each cited in `app.py` comments),
   and a 4-level FDHA logic-tree builder fed exclusively from the real model
   registry (`registry.json`), with live FDLT-001 weight validation and a
   single-file NRML preview. GUI-generated XML parses and validates cleanly
   with the engine's own `nrml_reader` / `validate_spec`.
2. **Run** — configuration summary, end-branch count, simulated progress,
   and a prominent DEMO MODE banner. Running is blocked while the logic
   tree is invalid.
3. **Results** — interactive hazard curves (mean + fractiles, log-log) and
   a displacement hazard map (mean + 5 fractile layers) with CSV downloads.

## Run locally

```bash
pip install -r webgui_demo/requirements.txt
streamlit run webgui_demo/app.py
```

The full oq-pfdha package is **not** required at runtime: model lists come
from the committed `registry.json` snapshot, results from the committed
CSVs.

### Regenerating the registry snapshot

If the model library changes, refresh the snapshot from the actual packages
(requires `pip install -e .` at the repository root):

```bash
python webgui_demo/generate_registry.py
```

### Regenerating the demo data

```bash
python openquake/fdha/test/benchmark/norcia_case3_iaea/run_test.py
# then re-copy the files listed in data/norcia_case3/PROVENANCE.md
```

## Deploy on Hugging Face Spaces (Docker)

1. Create a new Space → SDK: **Docker** (blank template).
2. Copy the contents of `webgui_demo/` to the Space repository root
   (`app.py`, `registry.json`, `data/`, `requirements.txt`, `Dockerfile`).
3. Push — Spaces builds the Dockerfile and serves the app on port 7860
   automatically.

Notes:
- The image is small (~400 MB): only Streamlit + plotting libraries, no
  OpenQuake stack.
- For the alternative *Streamlit SDK* Space (no Dockerfile), set
  `app_file: app.py` in the Space's README front matter and Spaces will use
  `requirements.txt` directly.

## License & citation

- License: **GNU AGPL v3.0 or later** (same as oq-pfdha; see `LICENSE` at
  the repository root).
- Please cite: Chen, Y.-S. (2025). *openquake.fdha: Python tools for
  probabilistic fault displacement hazard analysis* (v1.0.0) [Software].
  Istituto Nazionale di Oceanografia e di Geofisica Sperimentale (OGS).
  https://github.com/vup1120/oq-pfdha (see `CITATION.cff`).
- Benchmark reference: IAEA (2025), *Benchmarking Current Practices in
  Probabilistic Fault Displacement Hazard Analysis for Nuclear
  Installations*, IAEA-TECDOC-2092, DOI 10.61092/iaea.74us-dn4n.
