# [🟢 2d Bratu model](@id bratu)

```@contents
Pages = ["bratu.md"]
Depth = 3
```

We consider the problem of Mittelmann [^Farrell] [^Wouters]:

$$\Delta u + NL(\lambda,u) = 0$$

with Neumann boundary condition on $\Omega = (0,1)^2$ and where $NL(\lambda,u)\equiv-10(u-\lambda e^u)$. This is a good example to show how automatic branch switching works and also nonlinear deflation.

```@example BRATU
using CairoMakie
using Gridap
using Gridap.FESpaces

using GridapBifurcationKit
using BifurcationKit
const BK = BifurcationKit

BK.set_plot_backend!(BK.BK_Makie())
CairoMakie.activate!() # hide

# custom plot function to deal with Gridap
function plotgridap!(ax, x; k...)
    n = isqrt(length(x))
    heatmap!(ax, reshape(x, n, n); colormap = :viridis, k...)
end
plotgridap(x; k...) = (fig = Figure(); ax = Axis(fig[1,1]); plotgridap!(ax, x; k...); fig)
```

We are now ready to specify the problem using the setting of **Gridap.jl**: it allows to write the equations very closely to the mathematical formulation:

```@example BRATU
# discretisation
n = 40
domain = (0, 1, 0, 1)
cells = (n, n)
model = CartesianDiscreteModel(domain,cells)

# function spaces
order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, conformity=:H1,)#dirichlet_tags="boundary")
U = TrialFESpace(V)

Ω = Triangulation(model)
degree = 2*order
const dΩ = Measure(Ω, degree) # we make it const because it is used in res

# nonlinearity
NL(u) = exp(u)

# residual
res(u,p,v) = ∫( -∇(v)⋅∇(u) -  v ⋅ (u - p.λ ⋅ (NL ∘ u)) * 10 )*dΩ

# jacobian of the residual
jac(u,p,du,v) = ∫( -∇(v)⋅∇(du) - v ⋅ du ⋅ (1 - p.λ * ( NL ∘ u)) * 10 )*dΩ

# 3rd and 4th derivatives, used for aBS
d2res(u,p,du1,du2,v) = ∫( v ⋅ du1 ⋅ du2 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ
d3res(u,p,du1,du2,du3,v) = ∫( v ⋅ du1 ⋅ du2 ⋅ du3 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ

# example of initial guess
uh = zero(U)

# model parameter
par_bratu = (λ = 0.01,)

# weight for normbratu
const w = cumsum(ones(length(uh.free_values))) / length(uh.free_values)
w .= (1 .+ LinRange(-1,1,n+1)) * transpose(LinRange(-1,1,n+1)) |> vec
w .-= minimum(w)
normbratu(x) = norm(x .* w) / sqrt(length(x))

# problem definition
prob = GridapBifProblem(res, uh, par_bratu, V, U, dΩ, (@optic _.λ);
                jac = jac,
                # d2res = d2res,
                # d3res = d3res,
                plot_solution = (ax, x, p; ax1 = nothing, k...) -> plotgridap!(ax, x; k...),
                record_from_solution = (x, p;k...) -> normbratu(x))
```

In the same vein, we can continue this solution as function of $\lambda$:

```@example BRATU
optn = NewtonPar(eigsolver = EigArpack())
optc = ContinuationPar(p_max = 40., p_min = 0.01, ds = 0.01, max_steps = 1000, detect_bifurcation = 3, newton_options = optn, nev = 20, tol_stability = 1e-6, n_inversion = 6)
br = continuation(prob, PALC(tangent = Bordered()), optc;
	# plot = true,
	verbosity = 0,
	)
```

```@example BRATU
f,ax = BK.plot(br; dash_unstable_style = true)
f
```


## Automatic branch switching at simple branch points

We can compute the branch off the third bifurcation point:

```@example BRATU
br1 = continuation(br, 3,
        ContinuationPar(BK.getcontparams(br); ds = 0.005, dsmax = 0.05, max_steps = 140);
        # verbosity = 0, plot = true, 
        nev = 10,
        start_with_eigen = Val(false),
        scaleζ = norminf,
        callback_newton = BK.cbMaxNorm(10),
        )
```

You can also plot the two branches together:

```@example BRATU
f, ax = BK.plot(br,br1,plotfold=false; dash_unstable_style = true)
f
```

We continue our journey and compute the branch bifurcating of the first bifurcation point from the last branch we computed:

```@example BRATU
br2 = continuation(br1, 1,
        ContinuationPar(BK.getcontparams(br1);ds = 0.005, dsmax = 0.05, max_steps = 140);
        # verbosity = 0, plot = true, 
        nev = 10,
        start_with_eigen = Val(false),
        scaleζ = norminf,
        callback_newton = BK.cbMaxNorm(10),
        )
f, ax = BK.plot(br, br1, br2; dash_unstable_style = true)
f
```

## Automatic branch switching at the 2d-branch points

We now show how to perform automatic branch switching at the nonsimple branch points. However, we think it is important that the user is able to use the previous tools in case automatic branch switching fails. This is explained in the next sections.

The call for automatic branch switching is the same as in the case of simple branch points (see above) except that many branches are returned.

```@example BRATU
branches = continuation(br, 2,
        # verbosity = 0, plot = true,
        start_with_eigen = Val(false),
        usedeflation = true,
        verbosedeflation = false,
        callback_newton = BK.cbMaxNorm(10),
        )
```

You can plot the branches using

```@example BRATU
f, ax = BK.plot(br1, br2, branches..., br; dash_unstable_style = true)
f
```

## References

[^Farrell]:> Farrell, Patrick E., Casper H. L. Beentjes, and Ásgeir Birkisson. **The Computation of Disconnected Bifurcation Diagrams.** ArXiv:1603.00809 [Math], March 2, 2016.

[^Wouters]:> Wouters. **Automatic Exploration Techniques for the Numerical Continuation of Large–Scale Nonlinear Systems**, 2019.
