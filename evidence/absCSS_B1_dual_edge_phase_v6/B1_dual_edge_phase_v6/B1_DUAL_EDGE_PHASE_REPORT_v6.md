# B1 dual-edge phase-locked refinement — v6

## Purpose

Close the remaining grid-phase ambiguity in B1 by enforcing both compacton-edge phases exactly while h -> 0.  For prescribed (psi_L,psi_R),

    h_m = 2*pi / (m + psi_R - psi_L),
    L_m = N_m h_m,  N_m = round(40/h_m),

and x0 is chosen so that both edge phases are satisfied to floating-point roundoff.  All production diagnostics use T=1 and dt/h approximately 0.1.

The primary dynamic sequence uses (psi_L,psi_R)=(3/4,1/4), selected because previous B1 phase experiments showed strong late-time contamination near p=2.  A full 4x4 cardinal audit of the *initial* discrete defect is also included.  For p=2 the pair with the largest initial defect, (0,3/4), is independently propagated at two widely separated resolutions.

## M1, p=7/4, fixed (3/4,1/4)

| m | N | h | status | E_shape | A_rad(0.5) | dP/P | dH/H |
|---:|---:|---:|:---:|---:|---:|---:|---:|
| 63 | 398 | 0.10053096 | OK | 8.325190e-02 | 2.197299e-01 | 8.009080e-03 | -1.509289e-02 |
| 126 | 799 | 0.05006522 | OK | 2.196680e-02 | 7.629797e-02 | 5.563042e-04 | -1.043647e-03 |
| 188 | 1194 | 0.03351032 | OK | 9.419928e-03 | 3.482723e-02 | 1.088787e-04 | -2.042787e-04 |
| 251 | 1595 | 0.02508258 | OK | 4.790740e-03 | 1.748388e-02 | 2.334407e-05 | -4.382014e-05 |

Log-log fits over all four successful resolutions give

- shape error: E ~ h^2.0434;
- fixed-guard radiation amplitude: A_rad ~ h^1.7977.

The radiation amplitude decreases monotonically from 2.197299e-01 to 1.748388e-02.  No finite nonzero radiation plateau is visible.

## M1, p=2, fixed (3/4,1/4)

| m | N | h | status | E_shape | A_rad(0.5) | dP/P | dH/H |
|---:|---:|---:|:---:|---:|---:|---:|---:|
| 63 | 398 | 0.10053096 | FAIL | — | — | — | — |
| 126 | 799 | 0.05006522 | OK | 1.335465e-01 | 3.580246e-01 | 1.927815e-02 | -3.969300e-02 |
| 188 | 1194 | 0.03351032 | OK | 6.190675e-02 | 1.932203e-01 | 4.231873e-03 | -8.603247e-03 |
| 251 | 1595 | 0.02508258 | OK | 3.439021e-02 | 1.221421e-01 | 1.294266e-03 | -2.601900e-03 |

The coarse m=63/N=398 calculation fails at nonlinear step 92 (attempted t about 0.92).  The three finer levels all reach T=1.  Fits over those successful levels give

- shape error: E ~ h^1.9599;
- radiation amplitude: A_rad ~ h^1.5548.

Thus p=2 exhibits a genuine resolution threshold, but after the threshold both error and radiative tail decrease systematically.

## Independent p=2 phase pair: (0,3/4)

The cardinal initial-defect audit identifies (0,3/4) as the largest-F_inf pair for p=2 at every tested m.  Its dynamic spot-check is

| m | N | h | status | E_shape | A_rad(0.5) | dP/P | dH/H |
|---:|---:|---:|:---:|---:|---:|---:|---:|
| 126 | 807 | 0.04957148 | OK | 1.405856e-01 | 3.643067e-01 | 2.143100e-02 | -4.380109e-02 |
| 251 | 1603 | 0.02495803 | OK | 3.466076e-02 | 1.242465e-01 | 1.306306e-03 | -2.618084e-03 |

The two-level effective exponents are

- q_shape = 2.0405;
- q_Arad = 1.5676.

These agree closely with the primary phase-locked sequence and provide an independent check that its convergence is not a favorable-phase accident.

## M5 control, p=2, fixed (3/4,1/4)

M5 uses eta=0.05 h^3 -> 0.

| m | N | h | status | E_shape | A_rad(0.5) | dP/P | dH/H |
|---:|---:|---:|:---:|---:|---:|---:|---:|
| 63 | 398 | 0.10053096 | OK | 1.689735e-03 | 1.564346e-03 | -1.311236e-04 | 2.433135e-04 |
| 126 | 799 | 0.05006522 | OK | 6.061104e-04 | 8.379569e-04 | -3.156244e-05 | 6.079026e-05 |
| 188 | 1194 | 0.03351032 | OK | 3.348232e-04 | 4.905250e-04 | -1.412712e-05 | 2.755186e-05 |
| 251 | 1595 | 0.02508258 | OK | 2.177231e-04 | 3.949381e-04 | -7.897891e-06 | 1.550074e-05 |

Fits give

- shape error: E ~ h^1.4755;
- radiation amplitude: A_rad ~ h^1.0172.

M5 remains two to three orders of magnitude cleaner than M1 at the endpoint and has no nonlinear-solver failure in the sequence.

## Domain-length control

The exact two-edge construction lets L_m differ slightly from 40.  To test whether this creates the convergence trend, the m=126, p=2 calculation was repeated at the *same h and same edge phases* after adding one remote grid point, changing L from 40.002112035351 to 40.052177256922.

For M1 the relative changes were

- shape error: 8.657e-12;
- A_rad: 4.041e-12;
- drift P: 1.704e-11;
- drift H: 1.754e-11.

For M5 they are likewise at roundoff-level (1.463e-13 in shape error).  The O(h) adjustment of L is therefore not responsible for the observed convergence.

## Decision

**B1 = PASS_WITH_RESOLUTION_AND_PHASE_QUALIFICATION.**

1. M1 at p=7/4 converges cleanly when both edge phases are held fixed.
2. M1 at p=2 has a coarse-grid resolution threshold, but two independent fixed-phase geometries show approximately second-order shape convergence and about h^1.55--h^1.57 decay of the remote radiative amplitude once resolved.
3. M5 is the preferred production route: its radiation remains at 10^-3--10^-4 scale while eta -> 0.
4. The previous apparent nonconvergence from arbitrary N sequences was a mixture of true pre-asymptotic edge sensitivity and changing subcell geometry, not evidence for a nonzero limiting ripple.
5. This is a numerical benchmark conclusion, not a proof uniform over every continuous phase pair.

**B1 is frozen. B2 (exact compact dipole propagation) is unlocked.**
