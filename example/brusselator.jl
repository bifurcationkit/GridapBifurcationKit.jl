cd(@__DIR__)
using Pkg
pkg"activate ."

using Revise
using Plots, Gridap, SparseArrays#, Arpack
import BifurcationKit as BK
# using KrylovKit
using LinearAlgebra
using Gridap.FESpaces
using GridapBifurcationKit
#############################################
# discretisation
begin
Nx = 200
lx = 1
domain = (0, 1)
cells = (Nx)
model = CartesianDiscreteModel(domain, cells)

# function spaces
par_sh = (α = 2., β = 5.45, D1 = 0.008, D2 = 0.004, l = 0.3)

order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V0 = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags = "boundary")

Y = MultiFieldFESpace([V0, V0])
X = MultiFieldFESpace([TrialFESpace(V0, par_sh.α), TrialFESpace(V0, par_sh.β / par_sh.α)])

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω, degree)
end

function res((u,v),p,(V1, V2))
    NL1(X,Y) = X^2*Y - (p.β+1) * X + p.α
    NL2(X,Y) = p.β * X - X^2*Y

    ∫( -p.D1/p.l^2 * ∇(u)⋅∇(V1) + NL1∘(u,v) ⋅ V1 +
       -p.D2/p.l^2 * ∇(v)⋅∇(V2) + NL2∘(u,v) ⋅ V2 )*dΩ
end

function jac((u,v),p,(du, dv),(V1, V2))
    d1NL1(u,v) = 2*u*v - (p.β+1)
    d2NL1(u,v) = u*u
    d1NL2(u,v) = p.β - 2*u*v
    d2NL2(u,v) = -u*u

    ∫( -p.D1/p.l^2 * ∇(du)⋅∇(V1) + (d1NL1(u,v)*du + d2NL1(u,v)*dv) ⋅ V1 +
        -p.D2/p.l^2 * ∇(dv)⋅∇(V2) + (d1NL2(u,v)*du + d2NL2(u,v)*dv) ⋅ V2 )*dΩ
end

uh = zero(X)
uh.free_values .= 2
####################################################################################################
function plotsol!(X; nd = 60, k...)
    sol = X isa BK.BorderedArray ? X.u : X

    Plots.plot!(sol; k...)
end
plotsol2!(sol::Vector; k...) = (_uh = zero(X); _uh.free_values .= sol; plotsol!(_uh; k...))
plotsol2(sol; k...) = (plot();plotsol2!(sol; k...))
plotsol(sol; k...) = (plot();plotsol!(sol; k...))
####################################################################################################
using Statistics

recordSolSCH(x, p;k...) = (
            # u0 = reshape(x[1:length(x)÷2],Nx,Ny)[Nx÷2,Ny÷2],
            u0 = x[length(x)÷4],
            nrm =  norm(x[1:length(x)÷2], Inf),
            max = maximum(x[1:length(x)÷2]),
            # mean = mean(x[1:length(x)÷2]),
            zero = maximum(x[1:length(x)÷2]) - p,
            norminf = norm(x[1:length(x)÷2], Inf))

prob = GridapBifProblem(res, uh, par_sh, Y, X, dΩ, (BK.@optic _.l);
            jac = jac,
            plot_solution = (x,p;kp...) -> plotsol!(x; kp...),
            record_from_solution = recordSolSCH,
            # record_from_solution = (x, p) -> norm(x[1:length(x)÷2] .- p, 8),
            )

# Mass = GridapBifurcationKit.get_mass_matrix(prob)

eig = BK.EigArpack(0.1, :LM, tol = 1e-12)
# eig = BK.EigArnoldiMethod(sigma = 0.2, which = BK.LM(), tol = 1e-12)

sol = @time BK.solve(prob, BK.Newton(), BK.NewtonPar(verbose = true, tol=1e-11))

opts = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.01, p_max = 1.9, p_min= 0., detect_bifurcation = 3, nev = 30, plot_every_step = 20, newton_options = BK.NewtonPar(verbose = true, tol = 1e-11, eigsolver = eig), max_steps = 120, tol_stability = 1e-7, n_inversion = 4)
br = @time BK.continuation(prob, 
        # BK.PALC(),
        BK.AutoSwitch(),
        opts;
        plot = true,
        # verbosity = 2,
    )

plot(br)
plot(sol.u)

hopfpt = BK.get_normal_form(br, 1; verbose = true, scaleζ = BK.norminf, start_with_eigen = Val(false))
# ──▶ Hopf bifurcation point is: SuperCritical
# SuperCritical - Hopf bifurcation point at l ≈ 0.5130210968371925.
# Frequency ω ≈ 2.1395086289528544
# Period of the periodic orbit ≈ 2.936742213680523
# Normal form z⋅(iω + a⋅δl + b⋅|z|²):
# ┌─ a = 0.8770224100924121 + 0.566870399452485im
# └─ b = -0.4159343311347744 + 0.3798209689611216im
####################################################################################################
# index of the Hopf point in br.specialpoint
ind_hopf = 1

# newton iterations to compute the Hopf point
hopfpoint = BK.newton(br, ind_hopf; options = BK.NewtonPar(verbose=true, tol=1e-10), normN = BK.norminf, usehessian = false)
hopfpoint.u.p

optcdim2 = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, p_max = 6.5, p_min = 0.0, newton_options = opts.newton_options, detect_bifurcation = 0)

br_hopf = BK.continuation(br, ind_hopf, (BK.@optic _.β),
    optcdim2, verbosity = 2,
    # detection of codim 2 bifurcations with bisection
    detect_codim2_bifurcation = 0,
    # we update the Hopf problem at every continuation step
    update_minaug_every_step = 1,
    start_with_eigen = false,
    jacobian_ma = BK.MinAug(), # specific to large dimensions
    usehessian = false,
    plot = true,
    normC = BK.norminf)

scene = plot(br_hopf)
####################################################################################################
# automatic branch switching from Hopf point
opt_po = BK.NewtonPar(tol = 1e-8, verbose = true, max_iterations = 15)
opts_po_cont = BK.ContinuationPar(dsmin = 0.001,
        dsmax = 0.04, ds = 0.01,
        p_max = 2.0,
        max_steps = 100,
        newton_options = opt_po,
        plot_every_step = 1,
        detect_bifurcation = 3,
        nev = 11,
        tol_stability = 1e-6,
        )

hp = BK.get_normal_form(br,1; start_with_eigen = Val(false), scaleζ = BK.norminf)
# number of time slices for the periodic orbit
M = 41
probFD = BK.Trapeze(;M,
    # specific method for solving linear system
    # of Periodic orbits with trapeze method
    # You could use the default one :FullLU (slower here)
    jacobian = BK.FullSparseInplace(),
    N = length(sol.u),
    massmatrix = BK.getmassmatrix(prob, nothing, nothing),
    )

br_po = BK._continuation(
    # arguments for branch switching from the first
    # Hopf bifurcation point
    hp, BK.getprob(br),
    opts_po_cont, probFD;
    alg = BK.AutoSwitch(),
    verbosity = 3, plot = true,
    # arguments for continuation
    δp = 0.01,
    plot_solution = (x, p; kwargs...) -> begin
        _sol = BK.get_periodic_orbit(probFD, x, p.p)
        heatmap!(_sol.u'; ylabel="time", color=:viridis, kwargs...)
    end,
    normC = BK.norminf)

plot(br_po)