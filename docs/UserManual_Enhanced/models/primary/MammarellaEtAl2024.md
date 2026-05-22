## Mammarella et al. (2024) — Principal Surface Rupture Probability

Reference: Mammarella, L., et al. (2024). "Conditional probability of surface rupture: A numerical approach for principal faulting" Earthquake Spectra.

If `dip_mu` and/or `seismothickness` are not available in your source XML, add them explicitly in the INI configuration:

### Notes and cautions
- Invalid `MSR` or `HDD_str` values raise ValueError.
- Probability is clipped to [0, 1].
- Discretizations (N_LOGW, N_DIP, HDR grid, Z step) are internal implementation details and not configurable via the configuration file.
