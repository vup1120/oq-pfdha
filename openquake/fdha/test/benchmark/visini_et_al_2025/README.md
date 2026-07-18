This benchmark covers the Visini et al. (2025) Figure 13 and Figure 14
reproduction workflows, with source details recorded in
[REFERENCE.md](REFERENCE.md).

- **`fig13/`** reproduces the paper's Figure 13 decision-tree worked example
  (cases 1–3, conditional probabilities of exceedance) against the digitized
  published curves, with enforced tolerances in
  `fig13/test_fig13_reproduction.py`. The headline comparison figure is
  [Figures/visini2025_fig13_all_cases_ref_vs_impl.png](Figures/visini2025_fig13_all_cases_ref_vs_impl.png)
  (published curves dashed, oq-pfdha solid; median computed/reference ratios
  0.85–1.03, per-case stats in `fig13/fig13_agreement.json`). Each curve ends
  in a sharp roll-off to exactly zero - a genuine ±3σ truncation of the
  `Visini2025SecondaryFD` log-normal (matching the FDHLab reference
  convention), not a modelling error; see "Reading the truncation cliff" in
  `fig13/README.md`. Regenerate with:

  ```bash
  PYTHONPATH=. python openquake/fdha/test/benchmark/visini_et_al_2025/fig13/plot_fig13_cases_vs_reference.py
  ```

Run: `pytest openquake/fdha/test/benchmark/visini_et_al_2025 -q`
(add `-m slow` to include the Figure 13 reproduction cases).

Status: PASS for the Figure 13 reproduction benchmark. Figure 14 remains a
diagnostic workflow.
