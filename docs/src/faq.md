# FAQ

See also [FAQ of BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/faq/).

## Which mass matrix should I use for the stability of an incompressible flow?

For the incompressible Navier–Stokes equations, the pressure has no time derivative, so the
semi-discrete system is ``M\dot z = F(z, p)`` with a **singular** mass matrix whose pressure
block is zero. You should therefore pass the *velocity only* mass:

```julia
m0((u, p), (v, q)) = ∫( v ⊙ u ) * dΩ
prob = GridapBifProblem(res, u0, par, V, U, dΩ, lens; jac = jac, mass = m0)
```

Using the full L² mass (including ``\int p\,q``) or the identity matrix changes the spectrum
of the pencil and can spoil the detection of Hopf bifurcations. See
[Mass matrix](@ref mass-matrix).

## Which eigensolver should I use for a differential-algebraic problem?

`GridapBifProblem` is a `BifurcationKit.AbstractDAEBifProblem`: whatever eigensolver you put
in `NewtonPar(eigsolver = ...)` is automatically wrapped into `EigenDAE`, and the generalized
eigenvalue problem ``J\phi = \lambda M\phi`` is solved with the assembled mass matrix. You can
also build it explicitly with `EigenMassMatrix(M, eig)`. See [Eigen Solvers](@ref eigensolver-page).

## The eigenvalues returned by my shift-invert eigensolver look wrong

With a shift-invert eigensolver (`EigArpack`, `EigArnoldiMethod`), the option `which = LM()`
selects the eigenvalues **closest to the shift `sigma`**, not those with the largest real
part. To monitor the stability, choose a shift near the imaginary axis where the
bifurcation is expected (for example `sigma ≈ iω` for a Hopf point). See
[Eigen Solvers](@ref eigensolver-page).
