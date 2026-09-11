# [🟢 Brusselator (reaction–diffusion): Hopf bifurcation and periodic orbits](@id brusselator)

```@contents
Pages = ["brusselator.md"]
Depth = 3
```

We consider the Brusselator reaction–diffusion system on $\Omega=(0,1)$ with Dirichlet
boundary conditions:

```math
\begin{aligned}
\partial_t u &= \frac{D_1}{l^2}\,\partial_{xx} u + \alpha - (\beta+1)u + u^2 v,\\
\partial_t v &= \frac{D_2}{l^2}\,\partial_{xx} v + \beta u - u^2 v.
\end{aligned}
```

The homogeneous steady state $(u,v)=(\alpha,\beta/\alpha)$ is used as initial guess. We
continue it with respect to the parameter $l$, detect the Hopf bifurcation, compute its
normal form, continue the Hopf point in the second parameter $\beta$ and finally switch to
the branch of periodic orbits.

```@example BRUSSELATOR
using CairoMakie
CairoMakie.activate!()
using LinearAlgebra
using Gridap
using Gridap.FESpaces
using GridapBifurcationKit
using BifurcationKit
const BK = BifurcationKit
BK.set_plot_backend!(BK.BK_Makie())

# plot the two components of a solution (the first half of the dofs is u, the second is v)
function plotsol!(ax, x; k...)
    _x = x isa BK.BorderedArray ? x.u : x
    n = length(_x) ÷ 2
    lines!(ax, _x[1:n]; label = "u", k...)
    lines!(ax, _x[n+1:end]; label = "v", k...)
end
```

## FEM formulation

!!! tip "Discretisation"
    We use a small discretization in order to run this example on github CI...

We use a second order Lagrange interpolation and impose the initial guess on the boundary
through the `TrialFESpace`:

```@example BRUSSELATOR
Nx = 80
domain = (0, 1)
cells = (Nx)
model = CartesianDiscreteModel(domain, cells)

order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V0 = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags = "boundary")

Y = MultiFieldFESpace([V0, V0])

par = (α = 2., β = 5.45, D1 = 0.008, D2 = 0.004, l = 0.3)
X = MultiFieldFESpace([TrialFESpace(V0, par.α), TrialFESpace(V0, par.β / par.α)])

Ω = Triangulation(model)
degree = 2*order
const dΩ = Measure(Ω, degree)
```

## Weak forms

The residual and its jacobian are written close to the mathematical formulation:

```@example BRUSSELATOR
function res((u, v), p, (V1, V2))
    NL1(X, Y) = X^2 * Y - (p.β + 1) * X + p.α
    NL2(X, Y) = p.β * X - X^2 * Y

    ∫( -p.D1 / p.l^2 * ∇(u)⋅∇(V1) + NL1∘(u, v) ⋅ V1 +
       -p.D2 / p.l^2 * ∇(v)⋅∇(V2) + NL2∘(u, v) ⋅ V2 ) * dΩ
end

function jac((u, v), p, (du, dv), (V1, V2))
    d1NL1(u, v) = 2 * u * v - (p.β + 1)
    d2NL1(u, v) = u * u
    d1NL2(u, v) = p.β - 2 * u * v
    d2NL2(u, v) = -u * u

    ∫( -p.D1 / p.l^2 * ∇(du)⋅∇(V1) + (d1NL1(u, v)*du + d2NL1(u, v)*dv) ⋅ V1 +
        -p.D2 / p.l^2 * ∇(dv)⋅∇(V2) + (d1NL2(u, v)*du + d2NL2(u, v)*dv) ⋅ V2 ) * dΩ
end
```

## Bifurcation problem

```@example BRUSSELATOR
uh = zero(X)
uh.free_values .= 2

record_from_solution(x, p; k...) = (u0 = x[length(x) ÷ 4], nrm = norm(x[1:length(x) ÷ 2], Inf))

prob = GridapBifProblem(res, uh, par, Y, X, dΩ, (@optic _.l);
            jac = jac,
            plot_solution = (ax, x, p; ax1 = nothing, k...) -> plotsol!(ax, x; k...),
            record_from_solution = record_from_solution)
```

We first compute the steady state with a Newton solver:

```@example BRUSSELATOR
sol = BK.solve(prob, BK.Newton(), BK.NewtonPar(verbose = true, tol = 1e-11))
```

## Continuation and Hopf bifurcation

We now continue the steady state with respect to $l$, using the generalized eigensolver
associated to the mass matrix (handled automatically by `BifurcationKit` for an
`AbstractDAEBifProblem`):

```@example BRUSSELATOR
eig = BK.EigArpack(0.1, :LM, tol = 1e-12)

opts = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.01,
        p_max = 1.9, p_min = 0., detect_bifurcation = 3, nev = 30,
        newton_options = BK.NewtonPar(tol = 1e-11, eigsolver = eig),
        max_steps = 120, tol_stability = 1e-7, n_inversion = 4)
br = BK.continuation(prob, BK.AutoSwitch(), opts; verbosity = 0)

f, ax = BK.plot(br; dash_unstable_style = true)
f
```

The normal form of the first Hopf bifurcation point gives the frequency and the type of the
bifurcation:

```@example BRUSSELATOR
hopfpt = BK.get_normal_form(br, 1; verbose = false, scaleζ = BK.norminf, start_with_eigen = Val(false))
hopfpt.ω
```

## Hopf continuation (codim 2)

We continue the Hopf point in the second parameter $\beta$:

```@example BRUSSELATOR
optcdim2 = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds = 0.01,
        p_max = 6.5, p_min = 0., newton_options = opts.newton_options, detect_bifurcation = 0)

br_hopf = BK.continuation(br, 1, (@optic _.β), optcdim2;
        detect_codim2_bifurcation = 2,
        update_minaug_every_step = 1,
        start_with_eigen = false,
        jacobian_ma = BK.MinAug(),   # specific to large dimensions
        usehessian = false,
        normC = BK.norminf,
        verbosity = 0)

f, ax = BK.plot(br_hopf; dash_unstable_style = true)
f
```

## Periodic orbits

Finally, we switch to the branch of periodic orbits at the Hopf point. The mass matrix is
passed to the trapeze method to account for the differential-algebraic structure
$M \dot z = F(z,p)$:

```@example BRUSSELATOR
opt_po = BK.NewtonPar(tol = 1e-8, max_iterations = 15)
opts_po_cont = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.04, ds = 0.01,
        p_max = 2.0, max_steps = 50, newton_options = opt_po,
        detect_bifurcation = 3, nev = 11, tol_stability = 1e-4)

probFD = BK.Trapeze(; M = 41,
        jacobian = BK.FullSparseInplace(),
        massmatrix = BK.getmassmatrix(prob, nothing, nothing))

br_po = BK.continuation(br, 1, opts_po_cont, probFD;
        start_with_eigen = Val(false),
        δp = 0.01,
        plot_solution = (ax, x, p; ax1, iter, state, k...) -> begin
            _sol = BK.get_periodic_orbit(BK.getprob(iter), x, p.p)
            heatmap!(ax, _sol.u'; colormap = :viridis, k...)
        end,
        normC = BK.norminf)
```

```@example BRUSSELATOR
f, ax = BK.plot(br_po; dash_unstable_style = true)
f
```

```@example BRUSSELATOR
sol = BK.get_periodic_orbit(br_po, 20)
f = Figure(); ax = Axis(f[1,1], title = "Periodic solution", ylabel = "time")
heatmap!(ax, 1:2Nx, sol.t, sol.u)
f
```
