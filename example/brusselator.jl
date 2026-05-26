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
Nx = 400
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
    @assert false "WIP"    
    d1NL1(u,v) = -1 + 2*u*v + 2*p.σ*(u - 1/v)
    d2NL1(u,v) = u^2 + 2*p.σ*(u - 1/v)/v^2
    d1NL2(u,v) = -2*u*v - 2*p.σ*(u - 1/v)
    d2NL2(u,v) = -u^2 - 2*p.σ*(u - 1/v)/v^2

    ∫(  -∇(du)⋅∇(V1) - p.d * ∇(dv)⋅∇(V2) +
            d1NL1∘(u,v)⋅V1⋅du + d2NL1∘(u,v)⋅V1⋅dv +
            d1NL2∘(u,v)⋅V2⋅du + d2NL2∘(u,v)⋅V2⋅dv )*dΩ
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

prob = GridapBifProblem(res, uh, par_sh, Y, X, (BK.@optic _.l);
            # jac = jac,
            plot_solution = (x,p;kp...) -> plotsol!(x; kp...),
            record_from_solution = recordSolSCH,
            # record_from_solution = (x, p) -> norm(x[1:length(x)÷2] .- p, 8),
            )

Mass = GridapBifurcationKit.get_mass_matrix(prob, dΩ)

eig = BK.EigenMassMatrix(Mass, BK.EigArpack(0.1, :LM, tol = 1e-12))
# eig = BK.EigenMassMatrix(Mass, BK.EigArnoldiMethod(sigma = 0.2, which = BK.LM(), tol = 1e-12))

sol = @time BK.solve(prob, BK.Newton(), BK.NewtonPar(verbose = true, tol=1e-11))

opts = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.01, p_max = 1.9, p_min= 0., detect_bifurcation = 3, nev = 30, plot_every_step = 5, newton_options = BK.NewtonPar(verbose = true, tol = 1e-11, eigsolver = eig), max_steps = 120, tol_stability = 1e-7, n_inversion = 4)
br = @time BK.continuation(prob, BK.PALC(), opts;
        plot = true,
        verbosity = 2,
    )

plot(br)
plot(sol.u)

# c'est FAUX IL FAUT LEFT EV
hopfpt = BK.get_normal_form(br, 1; verbose = true, scaleζ = BK.norminf, autodiff = false)
# ──▶ Hopf bifurcation point is: SuperCritical
# SuperCritical - Hopf bifurcation point at l ≈ 0.513122950597178.
# Frequency ω ≈ 0.0026743076560430666
# Period of the periodic orbit ≈ 2349.4624087029138
# Normal form z⋅(iω + a⋅δl + b⋅|z|²):
# ┌─ a = 0.0010958820972717148 + 0.0007084183049488559im
# └─ b = -0.00038377805265247123 + 0.0005681844372716811im
####################################################################################################
# index of the Hopf point in br.specialpoint
ind_hopf = 1

# newton iterations to compute the Hopf point
hopfpoint = BK.newton(br, ind_hopf; options = BK.NewtonPar(verbose=true), normN = BK.norminf, usehessian = false)
hopfpoint.u.p

optcdim2 = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, p_max = 6.5, p_min = 0.0, newton_options = opts.newton_options, detect_bifurcation = 0)

br_hopf = BK.continuation(br, ind_hopf, (BK.@optic _.β),
    optcdim2, verbosity = 2,
    # detection of codim 2 bifurcations with bisection
    detect_codim2_bifurcation = 0,
    # we update the Hopf problem at every continuation step
    update_minaug_every_step = 1,
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
        max_steps = 30,
        newton_options = opt_po,
        plot_every_step = 1,
        detect_bifurcation = 0,
        nev = 11,
        tol_stability = 1e-6,
        )

# number of time slices for the periodic orbit
M = 41
probFD = BK.Trapeze(;M,
    # specific method for solving linear system
    # of Periodic orbits with trapeze method
    # You could use the default one :FullLU (slower here)
    jacobian = BK.FullSparseInplace(),
    N = length(sol.u),
    massmatrix = Mass
    )

br_po = BK.continuation(
    # arguments for branch switching from the first
    # Hopf bifurcation point
    br, 1,
    # arguments for continuation
    δp = 0.01, ampfactor = 1.91, use_normal_form = false,
    opts_po_cont, probFD;
    # regular options for continuation
    verbosity = 3, plot = true,
    autodiff_nf = false,
    plot_solution = (x, p; kwargs...) -> begin
        _sol = BK.get_periodic_orbit(probFD, x, p.p)
        heatmap!(_sol.u'; ylabel="time", color=:viridis, kwargs...)
    end,
    normC = BK.norminf)
