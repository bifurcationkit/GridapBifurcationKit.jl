# [Eigen solvers (Eig)](@id eigensolver-page)

See also [Eigen solvers](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/eigensolver/) in the `BifurcationKit.jl` documentation for more information, for example on how to implement your own eigensolver.

The eigensolver is specified through `NewtonPar(eigsolver = ...)` and is used to monitor the stability of the solutions during continuation (see [Bifurcation detection (codim 1)](@ref detection-page)).

## Generic eigen solvers

`BifurcationKit.jl` provides a list of eigensolvers which are directly usable here:

1. [`DefaultEig`](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/eigensolver/) for small problems (dense).
2. `EigArpack` for sparse problems.
3. `EigArnoldiMethod` for sparse problems.
4. `EigKrylovKit` for matrix-free problems.

## Generalized / DAE eigen solvers

For a differential-algebraic problem ``M\dot z = F(z, p)``, the stability is given by the
generalized eigenvalue problem

```math
J\,\phi = \lambda\, M\,\phi,
```

see [Mass matrix](@ref mass-matrix). When `GridapBifProblem` is used (it is an
`AbstractDAEBifProblem`), `BifurcationKit.jl` automatically wraps the user eigensolver into
`EigenDAE` and calls it with the mass matrix. One can also build the generalized solver
explicitly:

1. `EigenMassMatrix(M, eig)` solves ``J\phi = \lambda M\phi`` using the eigensolver `eig`.
   For example

```julia
using BifurcationKit
M0 = GridapBifurcationKit.get_mass_matrix(prob)
eig = EigenMassMatrix(M0, EigArnoldiMethod(sigma = 0.1, which = LM()))
λ, vp, conv, iter = eig(J, 10)
```

2. `EigenDAE(eig)` is the wrapper used internally for `AbstractDAEBifProblem`. Passing a
   plain eigensolver to `NewtonPar` is enough, the mass matrix is added automatically.

!!! tip "Choosing the shift"
    With a shift-invert eigensolver (`EigArpack`, `EigArnoldiMethod`), `which = LM()` returns
    the eigenvalues **closest to the shift `sigma`**, not those with the largest real part.
    For stability, choose a shift close to the imaginary axis where the bifurcation is
    expected (for example `sigma ≈ iω` for a Hopf bifurcation).

!!! warning "Singular mass matrix"
    When ``M`` is singular (*e.g.* a zero pressure block for incompressible flows), only the
    **finite** eigenvalues of the pencil are physical. The infinite eigenvalues associated
    with the algebraic constraints should be discarded (they are usually harmless when using
    a well chosen shift).
