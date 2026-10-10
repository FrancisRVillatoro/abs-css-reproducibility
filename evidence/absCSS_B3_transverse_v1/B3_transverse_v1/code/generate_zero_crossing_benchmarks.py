#!/usr/bin/env python3
"""Generate the initial-time transverse-zero regularisation benchmark.

For the Hamiltonian smoothing
    A_{eps,r}(u) = (u^(2r)+eps^(2r))^(p/(2r)),
and the imposed local crossing u=a*x, the regularised CSS dispersive-flux
derivative is

    K = p a^3 (a x)^(2r-2) ((a x)^(2r)+eps^(2r))^(p/(2r)-2)
        * [(2r-1) eps^(2r) + (p-1)(a x)^(2r)].

The script verifies the initial-time laws
    peak |K| ~ eps^(p-2),
    matching/core width ~ eps/|a|,
    integral over |x| <= Z*eps/|a| ~ eps^(p-1),
without making any claim about the subsequent PDE dynamics.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover
    raise SystemExit("matplotlib is required to generate the benchmark figures") from exc

DEFAULT_P = (1.0, 1.25, 1.5, 1.75)
DEFAULT_R = (1, 2)
DEFAULT_EPS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3)


def kernel(x: np.ndarray, p: float, eps: float, r: int, slope: float = 1.0) -> np.ndarray:
    if p <= 0 or eps <= 0 or r < 1 or slope == 0:
        raise ValueError("Require p>0, eps>0, integer r>=1, and nonzero slope")
    y = slope * x
    S = y ** (2 * r) + eps ** (2 * r)
    return (
        p
        * slope**3
        * y ** (2 * r - 2)
        * S ** (p / (2 * r) - 2)
        * ((2 * r - 1) * eps ** (2 * r) + (p - 1) * y ** (2 * r))
    )


def normalized_kernel(z: np.ndarray, p: float, r: int) -> np.ndarray:
    # K(eps*z/a)/(a^3 eps^(p-2)); independent of eps and positive slope a.
    S = z ** (2 * r) + 1.0
    return (
        p
        * z ** (2 * r - 2)
        * S ** (p / (2 * r) - 2)
        * ((2 * r - 1) + (p - 1) * z ** (2 * r))
    )


def fit_power(eps: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    coeff = np.polyfit(np.log(eps), np.log(np.abs(values)), 1)
    slope, intercept = float(coeff[0]), float(coeff[1])
    return slope, math.exp(intercept)


def peak_constant(p: float, r: int) -> float:
    z_linear = np.linspace(0.0, 20.0, 250_001)
    z_log = np.logspace(math.log10(20.0), 4.0, 50_000)
    z = np.concatenate((z_linear, z_log[1:]))
    phi = np.abs(normalized_kernel(z, p, r))
    return float(np.nanmax(phi))


def transition_constant(p: float, r: int, tolerance: float = 0.10) -> float:
    """One-sided z where the regularised kernel matches the outer singular law.

    For p>1 the outer law is p(p-1) z^(p-2).  We return the final z at
    which the relative error exceeds `tolerance`; physical width is z*eps/|a|.
    """
    if p <= 1:
        # Use an outer half-maximum width for the p=1 delta kernel.
        z = np.linspace(0.0, 50.0, 500_001)
        phi = np.abs(normalized_kernel(z, p, r))
        half = 0.5 * np.nanmax(phi)
        idx = np.flatnonzero(phi >= half)
        return float(z[idx[-1]]) if idx.size else float("nan")

    z = np.logspace(-5.0, 6.0, 300_000)
    phi = normalized_kernel(z, p, r)
    outer = p * (p - 1) * z ** (p - 2)
    rel = np.abs(phi - outer) / np.maximum(np.abs(outer), np.finfo(float).tiny)
    bad = np.flatnonzero(rel > tolerance)
    if not bad.size:
        return float(z[0])
    if bad[-1] == len(z) - 1:
        raise RuntimeError(f"z range insufficient for p={p}, r={r}")
    return float(z[bad[-1] + 1])


def scaled_window_integral(p: float, r: int, zmax: float = 4.0) -> float:
    z = np.linspace(-zmax, zmax, 400_001)
    return float(np.trapezoid(normalized_kernel(z, p, r), z))


def p1_total_mass(r: int) -> float:
    # A symmetric hybrid grid resolves the core and algebraic tails.
    zp = np.concatenate((np.linspace(0.0, 20.0, 300_001), np.logspace(math.log10(20.0), 6.0, 200_000)[1:]))
    phi = normalized_kernel(zp, 1.0, r)
    return float(2.0 * np.trapezoid(phi, zp))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("zero_crossing_benchmarks"))
    parser.add_argument("--p", nargs="*", type=float, default=list(DEFAULT_P))
    parser.add_argument("--r", nargs="*", type=int, default=list(DEFAULT_R))
    parser.add_argument("--eps", nargs="*", type=float, default=list(DEFAULT_EPS))
    parser.add_argument("--slope", type=float, default=1.0)
    parser.add_argument("--central-z", type=float, default=4.0)
    args = parser.parse_args()

    if any(p <= 0 for p in args.p):
        raise ValueError("All p values must be positive")
    if any(r < 1 for r in args.r):
        raise ValueError("All r values must be positive integers")
    if any(eps <= 0 for eps in args.eps):
        raise ValueError("All epsilon values must be positive")
    if args.slope <= 0:
        raise ValueError("This benchmark currently assumes a positive slope")

    outdir = args.outdir
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    eps_arr = np.array(sorted(args.eps, reverse=True), dtype=float)
    rows: list[dict[str, float | int]] = []
    fits: list[dict[str, float | int]] = []

    for r in args.r:
        for p in args.p:
            peak_c = peak_constant(p, r)
            width_c = transition_constant(p, r)
            mass_c = scaled_window_integral(p, r, args.central_z)

            peaks = []
            widths = []
            masses = []
            for eps in eps_arr:
                peak = abs(args.slope) ** 3 * eps ** (p - 2) * peak_c
                width = width_c * eps / abs(args.slope)
                mass = args.slope**2 * eps ** (p - 1) * mass_c
                peaks.append(peak)
                widths.append(width)
                masses.append(mass)
                rows.append(
                    {
                        "p": p,
                        "regularisation_r": r,
                        "epsilon": float(eps),
                        "slope_a": args.slope,
                        "peak_abs_K": peak,
                        "transition_or_halfmax_radius": width,
                        f"integral_absx_le_{args.central_z:g}_eps_over_a": mass,
                    }
                )

            peak_exp, peak_pref = fit_power(eps_arr, np.array(peaks))
            width_exp, width_pref = fit_power(eps_arr, np.array(widths))
            mass_exp, mass_pref = fit_power(eps_arr, np.array(masses))
            record: dict[str, float | int] = {
                "p": p,
                "regularisation_r": r,
                "fitted_peak_exponent": peak_exp,
                "predicted_peak_exponent": p - 2,
                "fitted_width_exponent": width_exp,
                "predicted_width_exponent": 1.0,
                "fitted_central_integral_exponent": mass_exp,
                "predicted_central_integral_exponent": p - 1,
                "peak_prefactor": peak_pref,
                "width_prefactor": width_pref,
                "central_integral_prefactor": mass_pref,
            }
            if abs(p - 1.0) < 1e-12:
                record["p1_total_normalized_mass"] = p1_total_mass(r)
                record["p1_predicted_total_normalized_mass"] = 2.0
            fits.append(record)

    csv_path = outdir / "initial_time_scaling.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fit_path = outdir / "fitted_exponents.json"
    fit_path.write_text(json.dumps(fits, indent=2), encoding="utf-8")

    # Figure 1: peak scaling.  Matplotlib default colours/styles are retained.
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for r in args.r:
        for p in args.p:
            subset = [row for row in rows if row["regularisation_r"] == r and row["p"] == p]
            xeps = np.array([row["epsilon"] for row in subset], dtype=float)
            ypeak = np.array([row["peak_abs_K"] for row in subset], dtype=float)
            ax.loglog(xeps, ypeak, marker="o", label=f"p={p:g}, r={r}")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$\max_x |D_x B_{p,\varepsilon,r}|$")
    ax.set_title("Initial-time transverse-zero peak scaling")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "zero_crossing_peak_scaling.png", dpi=180)
    plt.close(fig)

    # Figure 2: normalised kernels for representative p values.
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    zplot = np.linspace(-8.0, 8.0, 5001)
    representative = [p for p in (1.0, 1.5) if any(abs(p - x) < 1e-12 for x in args.p)]
    if not representative:
        representative = args.p[: min(2, len(args.p))]
    for r in args.r:
        for p in representative:
            ax.plot(zplot, normalized_kernel(zplot, p, r), label=f"p={p:g}, r={r}")
    ax.set_xlabel(r"$z=a x/\varepsilon$")
    ax.set_ylabel(r"$D_xB/(a^3\varepsilon^{p-2})$")
    ax.set_title("Collapse of normalised transverse-zero kernels")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figdir / "zero_crossing_normalized_kernels.png", dpi=180)
    plt.close(fig)

    # Figure 3: transition/half-maximum radius scaling.
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for r in args.r:
        for p in args.p:
            subset = [row for row in rows if row["regularisation_r"] == r and row["p"] == p]
            xeps = np.array([row["epsilon"] for row in subset], dtype=float)
            yw = np.array([row["transition_or_halfmax_radius"] for row in subset], dtype=float)
            ax.loglog(xeps, yw, marker="o", label=f"p={p:g}, r={r}")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("one-sided transition/core radius")
    ax.set_title("Initial-time transverse-zero width scaling")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "zero_crossing_width_scaling.png", dpi=180)
    plt.close(fig)

    manifest = {
        "status": "initial-time algebraic benchmark; not a PDE evolution",
        "p_values": args.p,
        "regularisation_r_values": args.r,
        "epsilon_values": list(eps_arr),
        "slope_a": args.slope,
        "central_window": f"|x| <= {args.central_z:g} epsilon/|a|",
        "predictions": {
            "peak_exponent": "p-2",
            "width_exponent": 1,
            "central_window_integral_exponent": "p-1",
            "p_equals_1_total_mass": "2 a |a|",
        },
        "files": {
            "data": csv_path.name,
            "fits": fit_path.name,
            "figures": [p.name for p in sorted(figdir.glob("*.png"))],
        },
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    max_peak_error = max(abs(float(item["fitted_peak_exponent"]) - float(item["predicted_peak_exponent"])) for item in fits)
    max_width_error = max(abs(float(item["fitted_width_exponent"]) - 1.0) for item in fits)
    max_mass_error = max(abs(float(item["fitted_central_integral_exponent"]) - float(item["predicted_central_integral_exponent"])) for item in fits)
    print(f"[OK] wrote {csv_path}")
    print(f"[OK] maximum exponent errors: peak={max_peak_error:.3e}, width={max_width_error:.3e}, central integral={max_mass_error:.3e}")
    for item in fits:
        if item["p"] == 1.0:
            print(f"[OK] p=1, r={item['regularisation_r']}: total normalised mass={item['p1_total_normalized_mass']:.12g}")


if __name__ == "__main__":
    main()
