# Changelog

## v1.0.0 — 2026-10-10 (public release candidate)

- Repacks audited S1 evidence as **code and data only**; explicitly excludes manuscript/supplement sources and rendered figures.
- Preserves the original B1–B7 programs and data byte-for-byte.
- Keeps the recomputed interval spectral evidence and B5 refinements from revision r1.
- Redirects figure generation to ignored `outputs/figures/` and removes LaTeX supplement generation.
- Adds GitHub and Zenodo software citation metadata, license, manifest, and WSL publishing instructions.
- Does not rerun or alter the scientific PDE and interval computations at packaging time.
