# Absolute CSS: reproducibility code and data

**Version:** 1.0.0 · **Release type:** software and research data · **Date:** 2026-10-10  
**Author:** Francisco R. Villatoro (Universidad de Málaga)

This repository contains **only source code, numerical data, finite spectral certificates,
scientific verification scripts and provenance metadata** for the research project
*Weak compactons, singular zeros and dynamical realisations in an absolute
Cooper–Shepard–Sodano equation*.

**The paper itself is intentionally NOT part of this repository or the Zenodo deposit.**
There is no manuscript PDF, LaTeX source, submitted article, editorial correspondence
or typeset supplementary paper. Plots can be regenerated locally into an ignored folder;
no plot PDF or PNG is committed to the public release.

## Contents

| Path | Contents |
| --- | --- |
| `repro/code/` | Reproducers and auditors for rigorous spectral verification and numerical postprocessing; B5 extension drivers, interval C backend |
| `repro/inputs/` | Exact rational parameter/frequency subdivision meshes and provenance |
| `repro/results/spectral/` | Recomputed low/direct/WKB finite evidence, primitive tests and coverage audits |
| `repro/results/b5/` | Additional CSS and RH spatial refinement calculations, including recorded unsuccessful attempts |
| `repro/results/analysis/` | Source-derived CSVs, verified diagnostics and figure-plotting data |
| `repro/results/verification_logs/` | Recorded verification logs (not a substitute for a new local rerun) |
| `evidence/` | The extracted **original B1–B7 numerical archives**, with their programs, raw fields, reports and CSVs. Original data bytes are preserved. |
| `EXECUTION_LEDGER.json` | Provenance of additional calculations; failed RH attempts explicitly distinguished from accepted runs |
| `ORIGINAL_ARCHIVE_SHA256.json` | Original B1–B7 archive SHA-256 digests |
| `SOURCE_PROVENANCE.json` | Exact provenance of the repackage and non-scientific output-path adjustments |
| `SHA256SUMS.txt` | SHA-256 manifest for the committed release files; check **before** regenerating outputs |

### Historical reports

The original B1–B7 evidence contains historical reports and preliminary interpretations.
Some were superseded by the final revisions. They are retained **verbatim for traceability**,
not as competing claims. Use the current finite numerical evidence and the recorded controls.

## Reproducing in WSL (Ubuntu or Debian)

The directed floating-point audit expects **Linux x86-64 LP64** with a C compiler that
respects IEC 60559 binary64 directed rounding and shared `libmpfr.so.6` (MPFR 4.2.2
was used to produce this version). **Run this in WSL Linux, not native Windows Python.**

Requirements for Debian/Ubuntu (installed via your operating system, when needed):
`gcc`, `libmpfr6`, `python3`, `python3-venv`, and `git`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r repro/requirements.txt
python scripts/audit_package.py
```

Full scientific verification using the stored PDE fields:

```bash
python repro/code/reproduce.py --workers 3
```

Faster mode: **verify stored spectral coverage, without recomputing the 46,185 interval rectangles**:

```bash
python repro/code/reproduce.py --workers 3 --skip-spectral
```

To repeat the additional **B5 PDE simulations** as well (significantly longer):

```bash
python repro/code/reproduce.py --workers 3 --rerun-pde
```

The default rebuilds the finite interval certificate, checks symbolic and numerical
identities, postprocesses the original and new PDE data, and regenerates plots in
`outputs/figures/`. **It does not generate any manuscript or typeset supplement.**
The original benchmark evidence in `evidence/` is never overwritten.

> **Important:** `reproduce.py` overwrites certain *tracked numerical result files* in
> `repro/results/` during regeneration; that is expected and may make `git status`
> show modifications. Run in a disposable copy if you want to preserve a clean checkout.
> `SHA256SUMS.txt` describes the deposited release *before* such regeneration.

## Rigorous numerical evidence

The finite proof includes 70 low-frequency exponent boxes, 46,185 direct
parameter–frequency rectangles, and 125 high-frequency WKB cone verifications.
The corrected high-frequency enclosure handles the entire angle interval,
including its limiting value at zero. The directed-arithmetic backend includes
binary64/MPFR primitives, checks, exact rational test inputs and high-precision
spot checks. The mathematical meaning of the certificate is explained in the
associated research article; this repository carries the **independent finite
evidence and executable audit**, not the article.

The B5 results include three spatial resolutions for viscosity 0.003,
three for 0.001, and four for 0.0003. Spatial differences are empirical
resolution diagnostics, **not rigorous PDE error bounds**, and the results
do not establish a regularisation-independent long-time scattering limit.

The additional B5 computations were executed in a separate Linux x86-64
working container, **not on the Picasso supercomputer**. That distinction is
recorded in `EXECUTION_LEDGER.json` and `repro/ENVIRONMENT.json`.

## Checking integrity

```bash
python repro/code/check_manifest.py
python scripts/audit_package.py
```

The first command checks the SHA-256 hashes of the exact released files.
The second verifies that this checkout is **code/data-only**, checks the
historical B1–B7 SHA-256 metadata and all submitted new-run ledger hashes,
and validates the Zenodo/GitHub metadata. Neither command needs the internet.

## Citation and release policy

Version `v1.0.0` is a **proposed initial public release tag**, and no DOI has
been assigned to this package. Once GitHub and Zenodo have archived the
release, amend `CITATION.cff` or future versions with the actual DOI.

See `CITATION.cff` and `.zenodo.json` for citation and deposition metadata.
The repository license is MIT; **review the copyright/permission information,
including rights over historical contributed code, before publishing**.

## Deposit

The recommended path is to push this tree to a new GitHub repository, enable
its Zenodo integration, then create a GitHub **release** for `v1.0.0`.
Zenodo can then archive the release automatically. Do **not** also upload the
manual Zenodo ZIP into a separate Zenodo record unless you deliberately want
two distinct records and DOIs. A self-contained manual ZIP is provided as an
alternative when GitHub–Zenodo integration is not used.

See [`PUBLICACION_WSL.md`](PUBLICACION_WSL.md) for copy-paste commands.
