# B6 — CSS chase compacton–anticompacton collision

## Scope
The verified protocol contained five mandatory benchmarks; B6 is the first post-gate chase-collision benchmark.  We use p=3/2, c_fast=3, c_slow=1, initial center separation 8, moving frame c0=1.  First support contact occurs at t_touch=(8-2pi)/2=0.8584073464 and free-flight center crossing at t=4.

## Regularization path
The collision cannot be followed robustly with the B1 single-pulse path eta=0.05 h^3.  Empirically the interaction layer requires eta~C h^2.  We therefore test eta=C h^2, which still vanishes as h->0.  C=1.2 is sufficient for several phases/resolutions through t=4; C=1.5 is more phase-robust through T=6.

## Main results
1. There is no complete annihilation at the center-crossing time.  At t=4 both a positive and a negative macroscopic extremum survive.
2. Along eta=1.2 h^2, m=95 and m=100 give nearly identical t=4 extrema (about 5.87 and -2.1) and signed masses, whereas the coarse m=65 is visibly more dissipative.
3. Long-time runs to T=6 with eta=1.5 h^2 are reproducible for m=75 and m=85: u_max=5.0395 and 5.0442; u_min=-1.9385 and -1.9399.  Thus both signs re-emerge after the main overlap.
4. The T=6 quantitative outcome is not yet independent of C.  At m=75, C=1.0,1.2,1.5 gives u_max=5.2429,5.1093,5.0395 and u_min=-2.2179,-2.0336,-1.9385.  Hence topology is robust, but amplitudes/phase shifts are not yet extrapolated to eta->0.
5. The CSS mass C=integral u dx is conserved to ~5e-12 in the long runs.  P and H drift because hyperviscosity is finite; their drift grows with C, as expected.

## Verdict
B6 = PASS_TOPOLOGICAL_SURVIVAL / OPEN_QUANTITATIVE_ETA_TO_ZERO_LONG_TIME.
The old shorthand “annihilation” is not supported by these resolved runs: the compacton and anticompacton lose amplitude and exchange/radiate strongly during overlap, but two macroscopic signs survive and re-emerge.  A publication claim about asymptotic outgoing amplitudes or phase shifts requires an additional long-time double limit h->0 at fixed eta followed by eta->0.
