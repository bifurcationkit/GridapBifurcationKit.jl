cd(@__DIR__)
using Pkg
pkg"activate ."

using Revise, Plots
using Gridap, Setfield
using Gridap.FESpaces
using BifurcationKit, GridapBifurcationKit
plotgridap!(x; k...) = (n=Int(sqrt(length(x)));heatmap!(reshape(x,n,n); color=:viridis, k...))
plotgridap(x; k...) =( plot();plotgridap!(x; k...))
#############################################
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
dΩ = Measure(Ω, degree)

NL(u) = exp(u)
res(u, p, v) = ∫( -∇(v)⋅∇(u) -  v ⋅ (u - p.λ ⋅ (NL ∘ u)) * 10 )*dΩ
jac(u, p, du, v) = ∫( -∇(v)⋅∇(du) - v ⋅ du ⋅ (1 - p.λ *( NL ∘ u)) * 10 )*dΩ
d2res(u, p, du1, du2, v) = ∫( v ⋅ du1 ⋅ du2 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ
d3res(u, p, du1, du2, du3, v) = ∫( v ⋅ du1 ⋅ du2 ⋅ du3 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ

uh = zero(U)
par_bratu = (λ = 0.01,)

# weight for normbratu
const w = cumsum(ones(length(uh.free_values))) / length(uh.free_values)
w .= (1 .+ LinRange(-1,1,n+1)) * transpose(LinRange(-1,1,n+1)) |> vec
w .-= minimum(w)
normbratu(x) = norm(x .* w) / sqrt(length(x))

prob = GridapBifProblem(res, uh, par_bratu, V, U, (@lens _.λ);
                jac = jac,
                # d2res = d2res,
                # d3res = d3res,
                plot_solution = (x,p; k...) -> plotgridap!(x;  k...),
                record_from_solution = (x, p) -> normbratu(x))

# factorize leads pivots issues, better use LU factorization here
optn = NewtonPar(eigsolver = EigArpack(), linsolver = DefaultLS(false))#EigKrylovKit(dim = 100))
sol = newton(prob, NewtonPar(optn; verbose = true))

opts = ContinuationPar(p_max = 40., p_min = 0.01, ds = 0.01, max_steps = 1000, detect_bifurcation = 3, newton_options = optn, nev = 20, tol_stability = 1e-6, n_inversion = 6, dsmin_bisection = 1e-17, max_bisection_steps = 25, tol_bisection_eigenvalue = 1e-19)
br = continuation(prob, PALC(tangent = Bordered()), opts;
    plot = true,
    verbosity = 0,
    )

plot(br)

nf = get_normal_form(br, 2; verbose = true, scaleζ = norminf)
####################################################################################################
br1 = continuation(br, 3,
        setproperties(opts; ds = 0.005, dsmax = 0.05, max_steps = 140, detect_bifurcation = 3);
        verbosity = 3, plot = true, nev = 10,
        usedeflation = true,
        # scaleζ = norminf,
        callback_newton = BifurcationKit.cbMaxNorm(100),
        )

plot(br, br1)

br2 = continuation(br1, 3,
        setproperties(opts;ds = 0.005, dsmax = 0.05, max_steps = 140, detect_bifurcation = 3);
        verbosity = 0, plot = true, nev = 10,
        usedeflation = true,
        callback_newton = BifurcationKit.cbMaxNorm(100),
        )

plot(br, br1, br2, legend=false)

br3 = continuation(br, 2,
        setproperties(opts; ds = 0.005, dsmax = 0.05, max_steps = 140, detect_bifurcation = 0);
        verbosity = 0, plot = true,
        usedeflation = true,
        verbosedeflation = false,
        callback_newton = BifurcationKit.cbMaxNorm(100),
        )

plot(br, br1, br2, br3..., legend=false)
plot(br, br3...)
