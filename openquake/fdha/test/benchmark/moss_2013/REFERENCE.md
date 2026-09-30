# Reference for Moss et al. (2013)

## Cited reference

Moss, R. E. S., Stanton, K. V., and Buelna, M. I. (2013). The impact of
material stiffness on the likelihood of fault rupture propagating to the
ground surface. *Seismological Research Letters*, 84(3), 485-488.
DOI: 10.1785/0220110109. Open-access author copy distributed under the SSA
Open Access Policy.

Secondary source (redraw only):
Moss, R., Thompson, S., Kuo, C.-H., Younesi, K., and Baumont, D. (2022).
*Reverse Fault PFDHA*. Report GIRS-2022-05 (revised 1/17/2024).
DOI: 10.34948/N3F595.

## Source of the numerical values

- Moss et al. (2013) Eq. 1 and the Figure 3 legend: the four logistic
  regressions (reverse / strike-slip, stiff / soft at Vs30 = 600 m/s);
  checked exactly in `unit/test_moss2013_psr.py`.
- Moss et al. (2013) Figure 3 (PDF page 4): the four plotted curves,
  digitised here.
- GIRS-2022-05 Figure 3.2 (report p. 7, PDF page 20): the redraw, digitised
  here; its Eqs 3.4-3.5 restate the reverse coefficients.

## What this case validates

`Moss2013PrimarySR` for both faulting styles and both site classes: the
plotted curves are the legend equations (a logistic refit recovers every
coefficient within 1.3 %) and the library reproduces them.

## Tolerance and its justification

Maximum absolute difference 0.02 and median 0.01 in probability. One pixel
of the paper's 150 ppi figure is about 0.004 in probability; axis
calibration adds at most 0.0025.

## How to reproduce

```bash
python digitize_moss2013_fig3.py /path/to/Mossetal.2013-OpenAccess.pdf
python digitize_girs_fig3_2.py /path/to/Moss_REV_1.17.2024.pdf
python plot_comparison.py
pytest openquake/fdha/test/benchmark/moss_2013 -q
```

## Status

PASS on 2026-09-28 against the original Figure 3 (all four curves). The
GIRS-2022-05 reverse soft-soil redraw is a strict xfail (report erratum).
