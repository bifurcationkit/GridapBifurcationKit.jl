# Bifurcation Problem

```@contents
Pages = ["problems.md"]
Depth = 3
```

We refer to the [docs](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/) of `BifurcationKit.jl` for an in-depth description of the bifurcation problems. Here, we only focus on the ones related to `Gridap`.

The main structure is [`GridapBifProblem`](@ref) which encodes a system of PDEs
discretized with `Gridap.jl`. The weak form is written with the usual `Gridap` API, for
example

```julia
res((u, v), p, (V1, V2)) = ∫( -p.D1 * ∇(u)⋅∇(V1) + NL1∘(u, v) ⋅ V1 +
                              -p.D2 * ∇(v)⋅∇(V2) + NL2∘(u, v) ⋅ V2 ) * dΩ
```

and the problem is instantiated with

```julia
prob = GridapBifProblem(res, u0, par, V, U, dΩ, (@optic _.D1); jac = jac, mass = m0)
```

See the [tutorials](@ref) for complete examples.

```@docs
GridapBifProblem
```

## [Mass matrix](@id mass-matrix)

Many evolution PDEs write, after semi-discretization in space,

```math
M\frac{\mathrm{d}z}{\mathrm{d}t} = F(z, p)
```

where ``M`` is a **mass matrix** which can be *singular*. A typical example is the
incompressible Navier–Stokes equations: the velocity has a time derivative but the pressure
does not, so that the mass matrix has a zero block on the pressure dofs.

The stability of a stationary solution is obtained from the **generalized eigenvalue
problem**

```math
J\,\phi = \lambda\, M\,\phi
```

and not from the spectrum of ``J`` alone (which contains spurious modes coming from the
saddle-point structure). For this reason, `GridapBifProblem` is a subtype of
`BifurcationKit.AbstractDAEBifProblem` and the mass matrix is used automatically by the DAE
eigensolvers, see [Eigen Solvers](@ref eigensolver-page).

The mass matrix is assembled with

```@docs
GridapBifurcationKit.get_mass_matrix
```

By default the L² mass ``\int u\cdot v`` is used. For incompressible flows, one passes the
*velocity only* mass, which produces the required zero block on the pressure:

```julia
m0((u, p), (v, q)) = ∫( v ⊙ u ) * dΩ
prob = GridapBifProblem(res, u0, par, V, U, dΩ, lens; jac = jac, mass = m0)
```

The following `BifurcationKit` functions are specialized so that the mass matrix is
correctly taken into account:

- `BifurcationKit.is_mass_matrix_constant(prob)` returns `true`;
- `BifurcationKit.getmassmatrix(prob, x, p)` returns [`get_mass_matrix`](@ref)`(prob)`.
