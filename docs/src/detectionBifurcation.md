# [Bifurcation detection (codim 1)](@id detection-page)

The bifurcations are detected during a call to `br = continuation(prob, alg, contParams::ContinuationPar; kwargs...)` by turning on the following flags:

- `contParams.detect_bifurcation = 2` (for eigenvalues based bifurcations)
- `contParams.detect_event = 2` (for other bifurcations like inclinations)

The bifurcation points are located by looking at the spectrum **e.g.** by monitoring the unstable eigenvalues. The eigenvalue ``\lambda`` is declared unstable if `real(λ) > contParams.tol_stability`. The located bifurcation points are then returned in `br.specialpoint`.

## Precise detection of bifurcation points using Bisection

Note that the bifurcation points detected when `detect_bifurcation = 2` can be rather *crude* localization of the true bifurcation points. Indeed, we only signal that, in between two continuation steps *which can be large*, a (several) bifurcation has been detected. Hence, we only have a rough idea of where the bifurcation is located, unless your `dsmax` is very small... This can be improved as follows.

If you choose `detect_bifurcation = 3`, a **bisection algorithm** is used to locate the bifurcation points more precisely. It means that we recursively track down the change in stability. Some options in [`ContinuationPar`](@ref) control this behavior:

- `n_inversion`: number of sign inversions in the bisection algorithm
- `max_bisection_steps` maximum number of bisection steps
- `tol_bisection_eigenvalue` tolerance on real part of eigenvalue to detect bifurcation points in the bisection steps

If this is still not enough, you can use a Newton solver to locate them very precisely. See [Fold / Hopf Continuation](@ref).

!!! tip "Bisection mode"
    During the bisection, the eigensolvers are called like `eig(J, nev; bisection = true)` in order to be able to adapt the solver precision.

## Large scale computations

The user must specify the number of eigenvalues to be computed (like `nev = 10`) in the parameters `::ContinuationPar` passed to `continuation`. Note that `nev` is automatically incremented whenever a bifurcation point is detected [^1]. Also, there is an option in `::ContinuationPar` to save (or not) the eigenvectors. This can be useful in memory limited environments (like on GPUs).

[^1]: In this case, the Krylov dimension is not increased because the eigensolver could be a direct solver. You might want to increase this dimension using the callbacks in [`continuation`](@ref).

## List of detected bifurcation points

|Bifurcation|index used|
|---|---|
| Fold | fold |
| Hopf | hopf |
| Bifurcation point (single eigenvalue stability change, Fold or branch point) | bp |

## Eigensolver

The user must provide an eigensolver by setting `NewtonPar.eigsolver` inside the
`::ContinuationPar` passed to `continuation`. See [`NewtonPar`](@ref) and
[`ContinuationPar`](@ref) for more information, and [Eigen Solvers](@ref eigensolver-page)
for the list of available solvers.

!!! tip "DAE / mass matrix"
    When the problem is a differential-algebraic problem (a `GridapBifProblem`, which has a
    mass matrix), `BifurcationKit.jl` wraps the eigensolver into `EigenDAE` and solves the
    generalized eigenvalue problem ``J\phi = \lambda M\phi``. This is required to correctly
    detect Hopf bifurcations, see [Mass matrix](@ref mass-matrix).

## Generic bifurcation

By this we mean a change in the dimension of the Jacobian kernel. The detection of Branch point is done by analysis of the spectrum of the Jacobian.

The detection is triggered by setting `detect_bifurcation > 1` in the parameter `::ContinuationPar` passed to `continuation`.

## Fold bifurcation

The detection of **Fold** point is done by monitoring the monotonicity of the parameter.

The detection is triggered by setting `detect_fold = true` in the parameter `::ContinuationPar` passed to `continuation`. When a **Fold** is detected on a branch `br`, a point is added to `br.foldpoint` allowing for later refinement using the function `newtonFold`.

## Hopf bifurcation

The detection of **Hopf** point is done by analysis of the spectrum of the Jacobian (or of the pencil ``(J, M)`` for a differential-algebraic problem).

The detection is triggered by setting `detect_bifurcation > 1` in the parameter `::ContinuationPar` passed to `continuation`. When a **Hopf point** is detected, a point is added to `br.specialpoint` allowing for later refinement using the function `newton_hopf`.
