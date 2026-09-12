cd(@__DIR__)
using Pkg
pkg"activate ."

using Revise
using CairoMakie, Gridap, SparseArrays#, Arpack
import BifurcationKit as BK
using LinearAlgebra
using Gridap.FESpaces
using GridapBifurcationKit

Makie.inline!(true)
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
function plotsol!(ax, X; nd = 60, k...)
    sol = X
    lines!(ax, sol; k...)
    ax
end
plotsol2!(ax, sol::Vector; k...) = (_uh = zero(X); _uh.free_values .= sol; plotsol!(ax, _uh; k...))
plotsol2(sol; k...) = (f=Figure();plotsol2!(Axis(f[1,1]), sol; k...))
plotsol(sol; k...) = (f=Figure();plotsol!(Axis(f[1,1])sol; k...))
####################################################################################################
using Statistics

recordSolSCH(x, p; k...) = (
            # u0 = reshape(x[1:length(x)÷2],Nx,Ny)[Nx÷2,Ny÷2],
            u0 = x[length(x)÷4],
            nrm =  norm(x[1:length(x)÷2], Inf),
            max = maximum(x[1:length(x)÷2]),
            # mean = mean(x[1:length(x)÷2]),
            zero = maximum(x[1:length(x)÷2]) - p,
            norminf = norm(x[1:length(x)÷2], Inf))

prob = GridapBifProblem(res, uh, par_sh, Y, X, dΩ, (BK.@optic _.l);
            jac = jac,
            plot_solution = (ax,x,p;ax1,iter,state,kp...) -> plotsol!(ax,x; kp...),
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

BK.plot(br)[1]
lines(sol.u)

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
    # plot = true,
    normC = BK.norminf)

BK.plot(br_hopf)[1]
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
    plot_solution = (ax, x, p; ax1, kwargs...) -> begin
        _sol = BK.get_periodic_orbit(probFD, x, p.p)
        ax.ylabel="time"
        heatmap!(ax, _sol.u; colormap=:viridis)
    end,
    normC = BK.norminf)

plot(br_po)

####################################################################################################
# Transient simulation with TransientFEOperator
# ∂ₜu = D1/l² Δu + u²v - (β+1)u + α
# ∂ₜv = D2/l² Δv + βu - u²v
# i.e. we solve   ∫(∂ₜu⋅V1 + ∂ₜv⋅V2)dΩ - res((u,v), p, (V1,V2)) = 0
#
# The Hopf bifurcation of the steady branch is at l ≈ 0.513 (see continuation above).
# For l > 0.513 the steady state is unstable and a limit cycle develops: e.g. l = 0.6.
function simulate(p = (; par_sh..., l = 0.6);
                  Δt = 0.01, tF = 60.0, θ = 0.5, amp = 1e-1,
                  plot_every = 50, plot = true)
    # trial/test spaces for this parameter set (Dirichlet data depend on p)
    Xp = MultiFieldFESpace([TrialFESpace(V0, p.α), TrialFESpace(V0, p.β / p.α)])

    # steady state at p, then small deterministic perturbation of the free dofs
    op_FE = FEOperator((x, v) -> res(x, p, v),
                       (x, dx, v) -> jac(x, p, dx, v), Xp, Y)
    nls0  = NLSolver(LUSolver(); method = :newton, show_trace = false)
    x0    = solve(FESolver(nls0), op_FE)
    x0.free_values .+= amp .* sin.(0.7362 .* (1:length(x0.free_values)))

    # transient residual and jacobians: res(t, u, v), jac = ∂res/∂u, jac_t = ∂res/∂(∂ₜu)
    res_t(t, (u, v), (V1, V2)) = ∫( ∂t(u) ⋅ V1 + ∂t(v) ⋅ V2 )dΩ - res((u, v), p, (V1, V2))
    jac_x(t, (u, v), (du, dv), (V1, V2)) = (-1) * jac((u, v), p, (du, dv), (V1, V2))
    jac_t(t, (u, v), (du, dv), (V1, V2)) = ∫( du ⋅ V1 + dv ⋅ V2 )dΩ

    op_t     = TransientFEOperator(res_t, (jac_x, jac_t), Xp, Y)
    solver_t = ThetaMethod(nls0, Δt, θ) # better for stiff Brusselator
    # solver_t =  RungeKutta(nls0, Δt, :EXRK_RungeKutta_4_4)
    sol_t    = solve(solver_t, op_t, 0.0, tF, x0)

    if plot
        fig = Figure()
        ax  = Axis(fig[1, 1], title = "solution")
        ylims!(ax, (0,4))
        axm = Axis(fig[1, 2], xlabel = "t", ylabel = "max")
        ts   = Float64[]
        maxs = Float64[]
        for (n, (tn, uhn)) in enumerate(sol_t)
            x = get_free_dof_values(uhn)
            push!(ts, tn)
            push!(maxs, maximum(x[1:length(x) ÷ 2]))
            n % plot_every == 0 || continue
            empty!(ax)
            plotsol!(ax, x[1:length(x) ÷ 2])
            ax.title = "t = $(round(tn; digits = 2))"
            empty!(axm)
            lines!(axm, ts, maxs)
            display(fig)
        end
    end
    return sol_t
end

sol_ev = simulate((; par_sh..., l = 1.8); Δt = 0.01, tF = 220.0, plot_every = 30)

####################################################################################################
