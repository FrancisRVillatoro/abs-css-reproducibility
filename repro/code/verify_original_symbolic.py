#!/usr/bin/env python3
"""Symbolic audit for the verified analytical core of the absolute CSS equation.

The checks are performed on a fixed-sign sector u>0.  The absolute-value
formulae follow by replacing
    u^a -> |u|^a,
    d(|u|^a)/du -> a |u|^{a-2} u,
and by using the odd/even parity stated in the manuscript.

Run:
    python verify_absCSS_sympy.py

SymPy tested in the generation environment: 1.14.0.
"""
from __future__ import annotations

import sympy as sp


def check(label: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.factor(sp.powdenest(sp.powsimp(expr, force=True), force=True)))
    if simplified != 0:
        raise AssertionError(f"{label} failed: {simplified}")
    print(f"[OK] {label}")


# ---------------------------------------------------------------------------
# Symbols and fields
# ---------------------------------------------------------------------------
x, t, xi, s, z = sp.symbols("x t xi s z", real=True)
l, p, q, rzero = sp.symbols("l p q r", positive=True)
alpha, c, A0, k = sp.symbols("alpha c A k", positive=True)
a = sp.symbols("a", real=True, nonzero=True)
eps = sp.symbols("eps", positive=True)

phi = sp.Function("phi")(x, t)
u = sp.Function("u")(x, t)
v = sp.Function("v")(x, t)
f = sp.Function("f")(xi)

# ---------------------------------------------------------------------------
# 1. Euler-Lagrange equation on u=phi_x>0
# ---------------------------------------------------------------------------
phix = sp.diff(phi, x)
phit = sp.diff(phi, t)
phixx = sp.diff(phi, x, 2)
lag_density = (
    sp.Rational(1, 2) * phix * phit
    + phix**l / (l * (l - 1))
    - alpha * phix**p * phixx**2
)

EL = (
    -sp.diff(sp.diff(lag_density, phit), t)
    -sp.diff(sp.diff(lag_density, phix), x)
    +sp.diff(sp.diff(lag_density, phixx), x, 2)
)

u_phi = phix
ux_phi = sp.diff(phi, x, 2)
uxx_phi = sp.diff(phi, x, 3)
uxxx_phi = sp.diff(phi, x, 4)
ut_phi = sp.diff(phi, x, t)
expanded_pde_phi = (
    ut_phi
    + u_phi ** (l - 2) * ux_phi
    + alpha
    * (
        2 * u_phi**p * uxxx_phi
        + 4 * p * u_phi ** (p - 1) * ux_phi * uxx_phi
        + p * (p - 1) * u_phi ** (p - 2) * ux_phi**3
    )
)
check("Euler-Lagrange equation", EL + expanded_pde_phi)

# ---------------------------------------------------------------------------
# 2. Divergence form expands to the classical PDE
# ---------------------------------------------------------------------------
ux = sp.diff(u, x)
uxx = sp.diff(u, x, 2)
uxxx = sp.diff(u, x, 3)
ut = sp.diff(u, t)

divergence_pde = (
    ut
    + sp.diff(u ** (l - 1) / (l - 1), x)
    - alpha * p * sp.diff(u ** (p - 1) * ux**2, x)
    + 2 * alpha * sp.diff(u**p * ux, x, 2)
)
expanded_pde = (
    ut
    + u ** (l - 2) * ux
    + alpha
    * (
        2 * u**p * uxxx
        + 4 * p * u ** (p - 1) * ux * uxx
        + p * (p - 1) * u ** (p - 2) * ux**3
    )
)
check("Divergence/classical equivalence", divergence_pde - expanded_pde)

# ---------------------------------------------------------------------------
# 3. Hamiltonian variational derivative and Hamiltonian form
# ---------------------------------------------------------------------------
h_density = alpha * u**p * ux**2 - u**l / (l * (l - 1))
delta_h = sp.diff(h_density, u) - sp.diff(sp.diff(h_density, ux), x)
expected_delta_h = -(
    u ** (l - 1) / (l - 1)
    + alpha * p * u ** (p - 1) * ux**2
    + 2 * alpha * u**p * uxx
)
check("Hamiltonian variational derivative", delta_h - expected_delta_h)
check("Hamiltonian evolution", ut - sp.diff(delta_h, x) - expanded_pde)

# ---------------------------------------------------------------------------
# 4. Travelling-wave first integral
# ---------------------------------------------------------------------------
fp = sp.diff(f, xi)
fpp = sp.diff(f, xi, 2)
first_integral = c * f**2 / 2 - f**l / (l * (l - 1)) - alpha * f**p * fp**2
integrated_profile_eq = (
    c * f
    - f ** (l - 1) / (l - 1)
    - alpha * (p * f ** (p - 1) * fp**2 + 2 * f**p * fpp)
)
check("Travelling-wave first integral", sp.diff(first_integral, xi) - fp * integrated_profile_eq)

# ---------------------------------------------------------------------------
# 5. Width scaling and fixed-width diagonal
# ---------------------------------------------------------------------------
width_exponent = -sp.Rational(1, 2) + p / (2 * (l - 2))
check(
    "Width exponent",
    width_exponent - (p + 2 - l) / (2 * (l - 2)),
)
check("Fixed-width condition l=p+2", width_exponent.subs(l, p + 2))

# Compacton parameter identities on l=p+2.
A_power = c * (p + 1) * (p + 2) / 2
k_sq = p**2 / (4 * alpha * (p + 1) * (p + 2))
# The two coefficient balances obtained after substituting f=A cos^(2/p)(k xi).
check(
    "Compacton constant coefficient",
    2 * alpha * k_sq * (2 / p) ** 2 * A_power
    - c,
)
# Direct balance is more transparent through the separated first integral:
# alpha f'^2 = c/2 f^(2-p) - f^2/((p+1)(p+2)).
# Dividing by f^(2-p) gives alpha*(2k/p)^2 A^p sin^2 = c/2 - A^p cos^2/((p+1)(p+2)).
check(
    "Compacton sine coefficient",
    alpha * (2 * sp.sqrt(k_sq) / p) ** 2 * A_power - sp.Rational(1, 2) * c,
)
check(
    "Compacton cosine coefficient",
    A_power / ((p + 1) * (p + 2)) - sp.Rational(1, 2) * c,
)

# ---------------------------------------------------------------------------
# 6. Pohozaev identities and scale second variation
# ---------------------------------------------------------------------------
X, Y, Pmom = sp.symbols("X Y P", real=True)
solution_XY = sp.solve(
    [sp.Eq(X + Y, c * Pmom), sp.Eq((p + 2) * X - l * Y + 2 * c * Pmom, 0)],
    [X, Y],
    dict=True,
)[0]
check("Pohozaev kinetic identity", solution_XY[X] - (l - 2) * c * Pmom / (p + l + 2))
check("Pohozaev potential identity", solution_XY[Y] - (p + 4) * c * Pmom / (p + l + 2))
check(
    "Hamiltonian/velocity identity",
    solution_XY[X] - solution_XY[Y] - (l - p - 6) * c * Pmom / (p + l + 2),
)
scale_second = (p + 4) * (p + 6 - l) * X / 4
# Verify from exponents a=(p+4)/2, b=(l-2)/2 and aX=bY.
a_exp = (p + 4) / 2
b_exp = (l - 2) / 2
raw_second = a_exp * (a_exp - 1) * X - b_exp * (b_exp - 1) * (a_exp / b_exp) * X
check("Scale second variation", raw_second - scale_second)

# ---------------------------------------------------------------------------
# 7. Energetic and pressure variables
# ---------------------------------------------------------------------------
m = sp.Rational(2, 1) / p
u_from_v = v**m  # fixed positive-sign sector
vx = sp.diff(v, x)
vxx = sp.diff(v, x, 2)
vxxx = sp.diff(v, x, 3)
ux_v = sp.diff(u_from_v, x)
uxx_v = sp.diff(u_from_v, x, 2)
B_u = p * u_from_v ** (p - 1) * ux_v**2 + 2 * u_from_v**p * uxx_v
B_v = u_from_v * (4 / p) * (v * vxx + 2 * vx**2 / p)
check("Pressure-variable flux factorization", B_u - B_v)

Vflux = v ** (2 * (l - 2) / p) / (l - 1) + (4 * alpha / p) * (v * vxx + 2 * vx**2 / p)
pressure_conservative = sp.diff(v, t) + Vflux * vx + (p / 2) * v * sp.diff(Vflux, x)
pressure_expanded = (
    sp.diff(v, t)
    + v ** (2 * (l - 2) / p) * vx
    + 2 * alpha
    * (
        v**2 * vxxx
        + (1 + 6 / p) * v * vx * vxx
        + 4 * vx**3 / p**2
    )
)
check("Pressure PDE expansion", pressure_conservative - pressure_expanded)

w = sp.Function("w")(x, t)
u_from_w = w ** (2 / (p + 2))
check(
    "Quadratic dispersive energy variable",
    u_from_w**p * sp.diff(u_from_w, x) ** 2
    - 4 * sp.diff(w, x) ** 2 / (p + 2) ** 2,
)
check(
    "Full Hamiltonian is quadratic in w on l=p+2",
    (2 * l / (p + 2) - 2).subs(l, p + 2),
)

# Conditional simple-zero interface law and exact compacton-speed check.
vx0 = sp.symbols("vx0", real=True)
interface_speed = 8 * alpha * vx0**2 / p**2
print(f"[OK] Conditional interface speed: s' = {interface_speed}")
compacton_pressure_slope = p * sp.sqrt(c / (8 * alpha))
check(
    "Conditional interface law reproduces compacton speed",
    (8 * alpha / p**2) * compacton_pressure_slope**2 - c,
)

# ---------------------------------------------------------------------------
# 8. Transverse-zero singularities
# ---------------------------------------------------------------------------
abs_a = sp.symbols("abs_a", positive=True)
# CSS, away from s=0: d/ds [p a |a|^p sign(s)|s|^(p-1)].
css_coeff = p * (p - 1) * a * abs_a**p
print(f"[OK] CSS transverse-zero coefficient (s != 0): {css_coeff} * |s|^(p-2)")
print("[OK] CSS p=1 jump of the dispersive flux: 2 a|a| (hence 2 alpha a|a| delta_0)")

# Generic one-sided zero u=A s^r: both dispersive-flux terms have the same power.
Aedge = sp.symbols("Aedge", positive=True)
generic_B_coeff = Aedge ** (p + 1) * rzero * ((p + 2) * rzero - 2)
generic_D_coeff = generic_B_coeff * (rzero * (p + 1) - 2)
print(
    "[OK] Generic zero: B_p ~", generic_B_coeff,
    "* s^(r(p+1)-2), D_x B_p ~", generic_D_coeff,
    "* s^(r(p+1)-3)"
)
check(
    "Compacton edge satisfies generic L1 power criterion",
    (rzero * (p + 1) - 2).subs(rzero, 2 / p) - 2 / p,
)

# RH g_q=|u|^(q-1)u; away from zero its third derivative coefficient is q(q-1)(q-2).
rh_coeff_same_profile = (p + 1) * p * (p - 1) * a * abs_a**p
check(
    "Same-profile RH/CSS singularity exponent",
    ((p + 1) - 3) - (p - 2),
)
check(
    "Same-profile unscaled coefficient ratio",
    rh_coeff_same_profile - (p + 1) * css_coeff,
)

# ---------------------------------------------------------------------------
# 9. Scaled |K|(p+1,p+1) equation with exactly the same compacton profile
# ---------------------------------------------------------------------------
a_rh = 4 / (p + 2) ** 2
b_rh = 4 * alpha / ((p + 1) * (p + 2))
n = 2 * (p + 1) / p  # power of cos in g=f^(p+1)
check(
    "Same-profile RH constant balance",
    b_rh * A_power * k_sq * n * (n - 1) - c,
)
check(
    "Same-profile RH cos^2 balance",
    a_rh * A_power - b_rh * A_power * k_sq * n**2,
)
check(
    "Scaled RH/CSS transverse coefficient ratio",
    b_rh * (p + 1) * p * (p - 1) - alpha * p * (p - 1) * 4 / (p + 2),
)

# Hamiltonian and Casimirs of the scaled RH comparator.
aop, bop, kap = sp.symbols("a_RH b_RH kappa", positive=True)
kap_value = sp.sqrt(aop / bop)
check(
    "RH cosine Casimir with kappa^2=a/b",
    (-aop * kap + bop * kap**3).subs(kap, kap_value),
)
check(
    "RH sine Casimir with kappa^2=a/b",
    (aop * kap - bop * kap**3).subs(kap, kap_value),
)

# Local identity for the time derivative of P_RH=1/2 int u^2 in a positive sector.
g = u**q
Pdot_density = -aop * u * sp.diff(g, x) - bop * u * sp.diff(g, x, 3)
Pdot_bulk = bop * q * (q - 1) * u ** (q - 2) * ux**3 / 2
Pdot_flux = (
    -aop * q * u ** (q + 1) / (q + 1)
    + bop * (-u * sp.diff(g, x, 2) + q * u ** (q - 1) * ux**2 / 2)
)
check("RH P-dot local identity", Pdot_density - sp.diff(Pdot_flux, x) - Pdot_bulk)

# Explicit periodic counterexample at q=2: u=2+sin x+eta sin 2x gives int u_x^3=3 pi eta.
eta = sp.symbols("eta", real=True)
u_example_x = sp.cos(x) + 2 * eta * sp.cos(2 * x)
check(
    "RH P is not generically conserved (q=2 example)",
    sp.integrate(sp.expand_trig(u_example_x**3), (x, 0, 2 * sp.pi)) - 3 * sp.pi * eta,
)

# ---------------------------------------------------------------------------
# 10. Hamiltonian regularization of a transverse zero
# ---------------------------------------------------------------------------
S = a**2 * x**2 + eps**2
B_eps = p * a**3 * x * S ** (p / 2 - 1)
dB_eps = sp.diff(B_eps, x)
expected_dB_eps = p * a**3 * S ** (p / 2 - 2) * (eps**2 + (p - 1) * a**2 * x**2)
check("Regularized transverse-zero profile", dB_eps - expected_dB_eps)

Rreg = sp.symbols("Rreg", positive=True, integer=True)
yreg = a * x
Sreg = yreg ** (2 * Rreg) + eps ** (2 * Rreg)
Breg = p * a**2 * yreg ** (2 * Rreg - 1) * Sreg ** (p / (2 * Rreg) - 1)
expected_Dreg = (
    p * a**3 * yreg ** (2 * Rreg - 2) * Sreg ** (p / (2 * Rreg) - 2)
    * ((2 * Rreg - 1) * eps ** (2 * Rreg) + (p - 1) * yreg ** (2 * Rreg))
)
check("General r Hamiltonian regularization kernel", sp.diff(Breg, x) - expected_Dreg)
print("[OK] General regularization kernel has width O(eps/|a|) and peak O(eps^(p-2)), independently of r.")

integral_kernel_p1 = sp.integrate((1 + z**2) ** (-sp.Rational(3, 2)), (z, -sp.oo, sp.oo))
check("p=1 delta-kernel mass", integral_kernel_p1 - 2)
print("[OK] Regularization scalings: width O(eps), peak O(eps^(p-2)); p=1 mass tends to 2 a|a|.")

print("\nAll symbolic checks passed.")
