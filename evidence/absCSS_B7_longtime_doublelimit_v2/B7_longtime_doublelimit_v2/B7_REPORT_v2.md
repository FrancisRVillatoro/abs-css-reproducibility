# B7 v2 — long-time CSS chase: phase-locked continuation threshold

## Purpose
Close the remaining long-time quantitative scattering question after B6 by taking the limits in the intended order: spatial refinement at fixed hyperviscosity and only then vanishing hyperviscosity. The collision is the p=3/2 chase with c_fast=3, c_slow=1. The comparison time is tau=2 after the free center crossing.

## Geometry correction
A naive m-sequence changes the relative grid phase of the two compactons because 8/h is not locked when h=2*pi/m. B7 therefore uses a phase-locked separation

S_m=(n_m+delta) h,

with fixed delta and parity, S_m=8+O(h). The four compacton edges and both centers retain the same subcell geometry as h changes. The final comparison time is T_m=S_m/2+2, hence T_m -> 6.

## Main new result
After phase locking, the long-time branch still exhibits a finite-looking continuation threshold in eta. The cheap screen at dt/h=0.2 gives the brackets

| m | h | lower value that fails | value that passes |
|---:|---:|---:|---:|
| 95  | 0.06614 | 0.003 | 0.004 |
| 105 | 0.05984 | 0.005 | 0.006 |
| 115 | 0.05464 | 0.005 | 0.006 |
| 125 | 0.05027 | 0.005 | 0.006 |
| 135 | 0.04654 | 0.005 | 0.006 |

The m=115, eta=0.006 passing point was repeated with dt/h=0.1 and also reaches the full target time. Thus the plateau near eta ~ 5e-3--6e-3 is not produced solely by the coarser temporal screen.

This is qualitatively different from B1/B5 early-time thresholds proportional to powers of h: on the phase-locked long-time branch the admissible eta does not show evidence of tending to zero over the resolved range.

## Quantitative eta sensitivity
The outgoing extrema remain strongly eta-dependent. For example on the same phase-locked family:

- m=95, eta=0.004: u_max=5.4091, u_min=-2.4380;
- m=95, eta=0.005: u_max=5.2116, u_min=-2.1539;
- m=105, eta=0.006: u_max=5.0766, u_min=-1.9668;
- m=115, eta=0.006 (dt/h=0.1): u_max=5.0322, u_min=-1.9050;
- m=125, eta=0.006: u_max=4.9994, u_min=-1.8592;
- m=135, eta=0.006: u_max=4.9788, u_min=-1.8313.

At eta=0.006 the spatial trend is regular for m=105--135, but eta=0.005 loses the long-time branch for m=105,115(screen),125 and 135. Consequently the scattering amplitudes cannot yet be extrapolated to eta=0 along this branch.

## Continuation predictor test
A trajectory-continuation experiment used the converged eta=0.01 history as a Newton predictor for smaller eta. It did not remove the severe stiffness/branch loss near eta=0.0075--0.005. Hence the obstacle is not explained simply by restarting Newton from a poor local predictor.

## Interpretation
The robust B6 conclusion survives: both signs survive and reemerge. What B7 does *not* support is a quantitative eta->0 long-time scattering map. Over the resolved phase-locked sequence, the branch that reaches tau=2 appears to require an absolute hyperviscosity of order 5e-3--6e-3 rather than eta(h)->0.

This is not a proof that the regularized PDE has a positive critical eta, because branch-following of the fully discrete implicit problem can still fail before the PDE does. It is, however, sufficient to reject a numerical claim of a certified vanishing-eta scattering map at T~6 with the present method.

## Status
- PASS: phase-locked construction removes the relative-edge phase ambiguity.
- PASS: topological survival/reemergence from B6 remains robust.
- PASS: fixed eta=0.006 has a regular spatial trend on fine phase-locked meshes.
- OPEN/NEGATIVE: no numerically certified eta->0 long-time scattering map.
- STOP RULE: do not spend further production compute trying arbitrary smaller eta with the same time-step/Newton formulation; a different continuation formulation would constitute a new numerical-method project rather than completion of the present benchmark campaign.
