# B4 — Dynamic selection between M1–M5

## Status

**Final benchmark status:** `PASS_CONDITIONAL_SELECTION / FAIL_GLOBAL_ALL_REALIZATIONS`.

B4 was executed with the prescribed order of limits: time refinement first, then spatial refinement at fixed regularization, followed by removal of hyperviscosity or Hamiltonian regularization. The central exponent is p=3/2; cross-p screens cover p=1,5/4,3/2,7/4,2.

## 1. Compatible compacton data

For p=3/2 and T=0.2, M1 is second-order in h in the resolved sequence. At fixed fine grid, Hamiltonian regularizations M3 (r=1) and M4 (r=2) approach both each other and M1 as epsilon decreases. Power-law fits give

- D34 ~ epsilon^1.165
- D(M3,M1) ~ epsilon^1.124
- D(M4,M1) ~ epsilon^1.236.

M5 approaches M1 as eta->0; the last three points give an empirical exponent about 0.53. Time halving changes the field by only 1.68e-8 in the central M5 control.

## 2. Exact compact dipole

The same result holds for the dipole: D34 tends to zero with epsilon and M5 tends to M1 with eta. Independent B2 data already show direct M1–M2 agreement for p=3/2 with D12 ~ h^1.971. Thus on the compatible pressure-regular dipole, pressure and u formulations select the same weak travelling state.

## 3. Prepared transverse zero

This is the decisive non-compatible test. For p=3/2, L=8 and T=0.01, M3 and M4 converge to one another as epsilon is removed. With epsilon=0.1,0.06,0.03,0.025, the fit is

D34 ~ epsilon^1.103,

and D34(epsilon=0.025)=2.557167e-04.

The h-first M5 sequence converges at fixed eta and, only after this spatial limit is controlled, approaches the same central Hamiltonian branch within the remaining finite-epsilon/finite-eta uncertainty. The raw unregularized M1 sequence is much slower spatially for this nonsmooth state and is not used as the reference for declaring different limits.

The M3–M4 agreement is also observed at p=1 and p=5/4; thus the two Hamiltonian regularizations do not show evidence for distinct macroscopic selections.

## 4. Pressure formulation M2

The M2 Jacobian was independently checked by directional finite differences. At p=2, where v=u algebraically, scaled M2 and M1 agree with relative L2 difference 2.094e-09. This rules out an algebraic implementation error.

For p=3/2 the transverse zero maps to a cusp v~sgn(x)|x|^(3/4). Time refinement at N=800 shows that dt/h=0.1 is inadequate, but dt/h<=0.05 converges rapidly. At fixed eta_v=1e-5, a block-scaled Newton solver resolves N=800,1200,1600; the N=1200 to N=1600 field difference is 4.206e-05. N=400 does not resolve this regularized layer.

However, a uniform eta_v->0 continuation could not be certified: changing fixed eta_v produces increasingly stiff/cusp-sensitive pressure problems and requires eta-dependent temporal/spatial resolution. Therefore B4 does **not** establish an alternative nonzero M2 limit; it establishes a domain/regularity restriction of the pressure realization.

## 5. Scientific conclusion

The data support a common selected evolution among the two Hamiltonian regularizations M3/M4 and the vanishing-hyperviscosity route M5, and recover M1 on compatible compacton/dipole data. M2 agrees when its pressure variable is regular (including the exact p=2 control and the pressure-regular dipole), but a generic transverse zero with p<2 is cusp-like in v and is not covered by the same numerical realization without an additional pressure regularization.

Hence the correct conclusion is not “five equivalent discretizations”. It is:

> **Selection is robust within the admissible realization class, while the pressure formulation has a genuinely smaller regularity domain at generic sign changes.**

At p=1 the transverse-zero contribution is distributional; unregularized M1/M5 are not classical grid realizations of that object. M3/M4 still produce mutually convergent smooth approximations. Therefore a global all-p, all-five universal-selection statement is not well-posed.

## 6. Gate decision

- `PASS_HAMILTONIAN_REGULARIZATION_INDEPENDENCE`
- `PASS_COMPATIBLE_PROFILE_SELECTION`
- `PASS_HAMILTONIAN_VS_VANISHING_HYPERVISCOSITY` (central transverse case, within extrapolation uncertainty)
- `PASS_M2_ALGEBRAIC_CONTROL_p2`
- `DOMAIN_RESTRICTION_M2_GENERIC_TRANSVERSE_p<2`
- `FAIL_GLOBAL_ALL_REALIZATIONS_UNIVERSALITY`

This is a scientific result rather than a numerical failure. The next benchmark is B5: CSS versus RH with exactly the same compact profile. Collisions remain blocked until B5 is completed; the realization-dependent selection statement should be fixed before any collision campaign is interpreted.
