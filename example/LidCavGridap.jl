
using Revise
using Gridap
using GridapMakie, GLMakie
Makie.inline!(true)
###############

function streamfunction(Uh)
    # espace pour la stream function ψ (H¹₀)
    reffeψ = ReferenceFE(lagrangian, Float64, order)
    Vψ  = TestFESpace(model, reffeψ; conformity=:H1, labels=labels, dirichlet_tags=["diri0","diri1"])
    Uψ  = TrialFESpace(Vψ) 
    uh.free_values .= Uh                     # Uh = état complet (vitesse+pression)
    u = uh[1]                                # champ vitesse FE (uh[1])
    ex, ey = VectorValue(1.0,0.0), VectorValue(0.0,1.0)
    a(ψ, v) = ∫( -∇(ψ)⋅∇(v) )*dΩ
    # ∫ curl(u) v = ∫ (u⋅ey)ₓ v − (u⋅ex)ᵧ v  ⇔  Green : ∫ (ux·vy − uy·vx)
    l(v) = ∫( (u⋅ex)*(∇(v)⋅ey) − (u⋅ey)*(∇(v)⋅ex) )*dΩ
    opψ = AffineFEOperator(a, l, Uψ, Vψ)
    return solve(LinearFESolver(LUSolver()), opψ)
end

function plotsol!(ax, Uh, N = 50; k...)
    Xr = LinRange(0,1,N)
    Yr = LinRange(0,1,N)
    uhc = [evaluate(Uh, Gridap.Point(x, y)) for x in Xr, y in Yr]
    contourf!(ax, Xr, Yr, norm.(uhc); colormap = :vik, k...)
    arrows2d!(ax, Xr, Yr, getindex.(uhc,1), getindex.(uhc,2);
            lengthscale = 0.05, color = :black)
    ax
end

function plotsol(Uh, N = 50; k...)
    f = Figure()
    ax = Axis(f[1,1])
    plotsol!(ax, Uh, N; k...)
    f
end

function plotsol!(ax, Uh::Vector, _sol, N = 50; grid_layout_perso = nothing, k...)
    Xr = Yr = LinRange(0,1,N)
    _sol.free_values .= Uh
    uhc = [evaluate(_sol[1], Gridap.Point(x, y)) for x in Xr, y in Yr]
    ax.title="||u||"
    contourf!(ax, Xr, Yr, norm.(uhc); colormap = :bwr, levels = 10)
    contour!(ax, Xr, Yr, norm.(uhc); color = :black, labels = true, labelsize = 13, levels = 10)
    arrows2d!(ax, Xr, Yr, getindex.(uhc,1), getindex.(uhc,2);
            lengthscale = 0.05, color = :black)

    if isnothing(grid_layout_perso)==false
        st = streamfunction(Uh)
        ax2 = Axis(grid_layout_perso[1,2], title = "stream")
        uhc = [evaluate(st, Gridap.Point(x, y)) for x in Xr, y in Yr]
        contourf!(ax2, Xr, Yr, uhc; colormap = :vik, levels = 10)
        contour!(ax2, Xr, Yr, uhc; color = :black, labels = true, labelsize = 13, levels = 10)
    end
end
# empty!(ax);plotsol!(ax, soln.u, uh, N = 50;);f

###############
begin
n = 80
domain = (0,1,0,1)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri1",[6,])
add_tag_from_tags!(labels,"diri0",[1,2,3,4,5,7,8])

D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["diri0","diri1"])

reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)

uD0 = VectorValue(0,0)
uD1 = VectorValue(1,0)
U = TrialFESpace(V,[uD0,uD1])
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

Re = 500.0
conv(u,∇u) = -(∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

a((u,p),(v,q)) = ∫( -∇(v)⊙∇(u)/Re + (∇⋅v)*p - q*(∇⋅u) )dΩ
c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ
res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)
jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v)
op = FEOperator(res,jac,X,Y)
# algop = Gridap.FESpaces.get_algebraic_operator(op)
# res0(u) = Gridap.FESpaces.residual(algop, u)
end
####################################################################################################
using LineSearches: BackTracking
begin
nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

sol_gp = solve(solver,op) #uh, ph
Nu = length(sol_gp[1].free_values)

plotsol(sol_gp[1])
end
####################################################################################################
# rest(t,(u,p),(v,q)) = ∫( ∂t(u)*v )dΩ + res((u,p),(v,q))
# jact(t,(u,p),(du,dp),(v,q)) = jac((u,p),(du,dp),(v,q))
# jac_t(t,(u,p),(du,dp),(v,q)) = ∫( du⋅v )dΩ

# Xt = TransientMultiFieldFESpace([V, Q])
# Yt = TransientMultiFieldFESpace([U, P])
# opt = TransientFEOperator(rest,(jact,jac_t),Xt,Yt)
# linear_solver = LUSolver()
# Δt = 0.05
# θ = 0.5
# ode_solver = ThetaMethod(linear_solver,Δt,θ)
# Ut = TransientTrialFESpace(V,(x,t)->0)
# u₀ = sol#interpolate_everywhere(0.0,Ut(0.))
# t₀ = 0.0
# T = 10.0
# uₕₜ = solve(ode_solver,opt,t₀,T,sol)
# for (uₕ,t) in uₕₜ
#     # pvd[t] = createvtk(Ω,"poisson_transient_solution_$t"*".vtu",cellfields=["u"=>uₕ])
#     plotsol(uₕ)
# end
####################################################################################################
# bifurcation diagram
using GLMakie, BifurcationKit, GridapBifurcationKit
const BK = BifurcationKit
BK.set_plot_backend!(BK.BK_Makie())

begin
    conv_BK(u,∇u) = -(∇u')⋅u
    dconv_BK(du,∇du,u,∇u) = conv_BK(u,∇du)+conv_BK(du,∇u)

    a_BK((u,p),par,(v,q)) = ∫( -∇(v)⊙∇(u)/par.Re + (∇⋅v)*p - q*(∇⋅u) )dΩ
    c_BK(u,v) = ∫( v⊙(conv_BK∘(u,∇(u))) )dΩ
    dc_BK(u,du,v) = ∫( v⊙(dconv_BK∘(du,∇(du),u,∇(u))) )dΩ
    res_BK((u,p),par,(v,q)) = a_BK((u,p),par,(v,q)) + c_BK(u,v)
    jac_BK((u,p),par,(du,dp),(v,q)) = a_BK((du,dp),par,(v,q)) + dc_BK(u,du,v)
    m((u,p),(v,q)) = ∫( v⊙u + p*q )dΩ
    m0((u,p),(v,q)) = ∫( v⊙u  )dΩ
end

function record(Uh,p;k...)
    uh.free_values .= Uh
    u = uh[1]
    (;E = sum(∫( u ⊙ u  )dΩ)/2,)
end
# record(soln.u,1)

const uh = zero(Y)
uh.free_values .= rand(length(uh.free_values))
uh.free_values .= sol_gp.free_values
par_lid = (Re = 1000.0, a=1)
prob = GridapBifProblem(res_BK, uh, par_lid, Y, X, dΩ, (@optic _.Re); 
                        jac = jac_BK,
                        mass = m0,
                        record_from_solution = record,
                        plot_solution = (ax,x,p; ax1=nothing, k...) -> plotsol!(ax,x,uh,50; k...)
                        )


eig = EigArpack(sigma = 1., which = :LM) # you must chose a large enough sigma so that the imaginary part distance is not much
eig = EigArnoldiMethod(which = :LM, sigma = 1.)
optn = NewtonPar(eigsolver = eig)
soln = @time BifurcationKit.solve(prob, Newton(), NewtonPar(optn; verbose = true, tol = 1e-10))

# 10.1017/9781108863148 gives Hopf1 at Re = 8023, Hopf2 at Re = 8700
# VP1 = ±2.7762im
opts = ContinuationPar(p_max = 10000., p_min = 0.01, ds = 1., dsmax = 150., max_steps = 3000, detect_bifurcation = 0, newton_options = optn, nev = 300, tol_stability = 1e-6, n_inversion = 6, plot_every_step = 10, detect_event = 2, save_eigenvectors = false)
@reset opts.newton_options.verbose = true
br = @time continuation(prob, 
        # BK.AutoSwitch(tol_param = 0.1), 
        BK.Natural(),
        opts;
        plot = true, 
        verbosity = 3,
        event = BK.SaveAtEvent((7750, 8000, 8375, 8875, 9375, 9875,10375))
    )

BK.plot(br)[1]

hp = @time BK.get_normal_form(br, 1; start_with_eigen = Val(false),
					bls = BorderingBLS(optn.linsolver),
					)

_x = soln.u; _par = par_lid
_ind = 2;_x = br.specialpoint[_ind].x;_par = @set par_lid.Re = br.specialpoint[_ind].param
_J = @time BifurcationKit.jacobian(prob, _x, _par)
_M = BK.getmassmatrix(prob, _x, _par)
λs, Xv, = BK.gev(eig, _J, _M, 200); λs

# empty!(ax);plotsol!(ax, _x, uh, N = 50;);f
norm(_J * Xv - _M * Xv * Diagonal(λs), Inf)
####################################################################################################
# time solution
# Re = 8000.0
function simulate_direct(;Re = 7740., Δt = 0.01, θ = 1.0, t0 = 0., tF = 2.0, plot_every = 1)
    conv(u,∇u) = -(∇u')⋅u
    dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)
    a((u,p),(v,q)) = ∫( -∇(v)⊙∇(u)/Re + (∇⋅v)*p - q*(∇⋅u) )dΩ
    c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
    dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ
    res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)
    jac((u,p),(du,dp),(v,q)) = a((du,dp),(v,q)) + dc(u,du,v)
    mass_t(t, (dtu, dtp), (v, q)) = ∫( dtu⊙v )dΩ
    res_t(t, (u,p), (v,q)) = ∫( ∂t(u) ⊙ v )dΩ - res((u,p),(v,q))
    res_t_semi(t, (u,p), (v,q)) =  (-1) * res((u,p),(v,q))
    jac_a(t, (u,p), (du,dp), (v,q)) = (-1) * jac((u,p),(du,dp),(v,q))
    jac_t(t, (u,p), (du,dp), (v,q)) = mass_t(t, (du,dp), (v,q))
    E(u) = sum(∫( u ⊙ u  )dΩ)/2
    op_transient = TransientSemilinearFEOperator(mass_t, res_t_semi, jac_a, jac_t, X,Y;
                                                 constant_mass = true)
    # op_transient = TransientFEOperator(res_t, (jac_a, jac_t), X, Y)

    nls = NLSolver(show_trace=true,
                    method=:newton,
                    ftol = 1e-8,    # 1e-10 était trop serré pour le NS
                    xtol = 1e-10,   # inchangé
                    iterations = 30,  # 5 était trop bas
                    )
    solver_t = ThetaMethod(nls, Δt, θ)
    uh0 = sol_gp
    uh0.free_values .= hp.x0
    sol_t = solve(solver_t, op_transient, t0, tF, uh0)

    fig = Figure()
    ax  = Axis(fig[1, 1], title = "solution")
    axm = Axis(fig[1, 2], xlabel = "t")
    ts   = Float64[]
    data = Float64[]

    for (n, (tn, uhn)) in enumerate(sol_t)
        @info n tn
        push!(ts, tn)
        push!(data, E(uhn[1]))
        n % plot_every == 0 || continue
        empty!(ax)
        plotsol!(ax, uhn[1])
        ax.title = "t = $(round(tn; digits = 2)) - [$n]"
        empty!(axm)
        lines!(axm, ts, data)
        display(fig)
    end
end
simulate_direct(;Re = 7740., plot_every = 10, tF = 300)
