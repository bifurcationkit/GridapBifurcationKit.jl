# Detection of bifurcation points

The bifurcations are detected during a call to `br = continuation(prob, alg, contParams::ContinuationPar;kwargs...)` by turning on the following flags:

- `contParams.detect_bifurcation = 2` (for eigenvalues based bifurcations)
- `contParams.detect_event = 2` (for other bifurcations like inclinations)
    
## Precise detection of bifurcation points using Bisection    

Note that the bifurcation points detected when `detect_bifurcation = 2` can be rather *crude*  localization of the true bifurcation points. Indeed, we only signal that, in between two continuation steps *which can be large*, a (several) bifurcation has been detected. Hence, we only have a rough idea of where the bifurcation is located, unless your `dsmax` is very small... This can be improved as follows.

If you choose `detect_bifurcation = 3`, a **bisection algorithm** is used to locate the bifurcation points more precisely. It means that we recursively track down the change in stability. Some options in [`ContinuationPar`](@ref) control this behavior:

- `n_inversion`: number of sign inversions in the bisection algorithm
- `max_bisection_steps` maximum number of bisection steps
- `tol_bisection_eigenvalue` tolerance on real part of eigenvalue to detect bifurcation points in the bisection steps

!!! tip "Bisection mode"
    During the bisection, the eigensolvers are called like `eil(J, nev; bisection = true)` in order to be able to adapt the solver precision.
