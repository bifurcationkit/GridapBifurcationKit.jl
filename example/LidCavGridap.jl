
using Gridap
using GridapMakie, GLMakie
Makie.inline!(true)
###############
function plotsol(Uh, N = 30; k...)
	X = LinRange(0,1,N)
	Y = LinRange(0,1,N)
	uhc = [evaluate(Uh, Gridap.Point(x, y)) for x in X, y in Y]
	contourf(X,Y,norm.(uhc))
end

function plotsol!(ax, Uh::Vector, _sol, N = 30; k...)
	X = LinRange(0,1,N)
	Y = LinRange(0,1,N)
	_sol.free_values .= Uh
	uhc = [evaluate(_sol[1], Gridap.Point(x, y)) for x in X, y in Y]
	contourf!(ax, X,Y,norm.(uhc))
end

###############
begin
n = 90
domain = (0,1,0,1)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

for I in eachindex(model.grid.node_coords)
	@info I; model.grid.node_coords[I]
end

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
end

const Re = 100.0
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
####################################################################################################
using LineSearches: BackTracking
nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

sol = solve(solver,op) #uh, ph
Nu = length(sol[1].free_values)

plotsol(sol[1])

GLMakie.contour(Ωₕ, sol[1])
####################################################################################################
rest(t,(u,p),(v,q)) = ∫( ∂t(u)*v )dΩ + res((u,p),(v,q))
jact(t,(u,p),(du,dp),(v,q)) = jac((u,p),(du,dp),(v,q))
jac_t(t,(u,p),(du,dp),(v,q)) = ∫( du⋅v )dΩ

Xt = TransientMultiFieldFESpace([V, Q])
Yt = TransientMultiFieldFESpace([U, P])
opt = TransientFEOperator(rest,(jact,jac_t),Xt,Yt)
linear_solver = LUSolver()
Δt = 0.05
θ = 0.5
ode_solver = ThetaMethod(linear_solver,Δt,θ)
Ut = TransientTrialFESpace(V,(x,t)->0)
u₀ = sol#interpolate_everywhere(0.0,Ut(0.))
t₀ = 0.0
T = 10.0
uₕₜ = solve(ode_solver,opt,t₀,T,sol)
for (uₕ,t) in uₕₜ
    # pvd[t] = createvtk(Ω,"poisson_transient_solution_$t"*".vtu",cellfields=["u"=>uₕ])
	plotsol(uₕ)
end
####################################################################################################
# bifurcation diagram
using Plots, BifurcationKit, GridapBifurcationKit

convBK(u,∇u) = -(∇u')⋅u
dconvBK(du,∇du,u,∇u) = convBK(u,∇du)+convBK(du,∇u)

aBK((u,p),par,(v,q)) = ∫( -∇(v)⊙∇(u)/par.Re + (∇⋅v)*p - q*(∇⋅u) )dΩ
cBK(u,v) = ∫( v⊙(convBK∘(u,∇(u))) )dΩ
dcBK(u,du,v) = ∫( v⊙(dconvBK∘(du,∇(du),u,∇(u))) )dΩ
resBK((u,p),par,(v,q)) = aBK((u,p),par,(v,q)) + cBK(u,v)
jacBK((u,p),par,(du,dp),(v,q)) = aBK((du,dp),par,(v,q)) + dcBK(u,du,v)
m((u,p),(v,q)) = ∫( v⊙u + p*q )dΩ
m0((u,p),(v,q)) = ∫( v⊙u  )dΩ

uh = zero(Y)
uh.free_values .= rand(length(uh.free_values))
# uh.free_values .= sol.free_values
par_lid = (Re = 100.0, a=1)
prob = GridapBifProblem(resBK, uh, par_lid, Y, X, (@optic _.Re); 
						jac = jacBK,
						plot_solution = (ax,x,p; ax1=nothing, k...) -> plotsol!(ax,x,uh,30; k...))
# mass matrix
M = assemble_matrix(m,X,Y)
M0 = assemble_matrix(m0,X,Y)


optn = NewtonPar(eigsolver = EigArpack())
optn = NewtonPar(eigsolver = EigArnoldiMethod(sigma = 0.1))
# @set! optn.eigsolver = EigKrylovKit(dim = 200, verbose = 2)
# @set! optn.linsolver = GMRESKrylovKit(dim = 100, Pl = (UpperTriangular(_J)), verbose = 2)
soln = @time BifurcationKit.solve(prob, Newton(), NewtonPar(optn; verbose = true, tol = 1e-10))

norm( soln.u, Inf)
norm(sol.free_values - soln.u, Inf)

# sol.free_values .= soln

opts = ContinuationPar(p_max = 17000., p_min = 0.01, ds = 1., dsmax = 150., max_steps = 1000, detect_bifurcation = 0, newton_options = optn, nev = 20, tol_stability = 1e-6, n_inversion = 6, plot_every_step = 1, save_sol_every_step = 1)
@reset opts.newton_options.verbose = true
br = @time continuation(prob, PALC(), opts;
	plot = true, verbosity = 3,
	)


_J = BifurcationKit.jacobian(prob, soln.u, par_lid)
ind = 70
_J = prob(Val(:Jac), br.sol[ind].x, @set par_lid.Re = br.sol[ind].p)
println("\n\n--> parameter = ", br.sol[ind].p)
# _eigs = optn.eigsolver(-_J, 10)

_λs, _vps = @time optn.eigsolver(_J,10)
_λs, _vps = @time BifurcationKit.gev(optn.eigsolver,_J,M0,10)


br.sol[30].p


_J = prob(Val(:Jac), soln, par_lid)

_J+_J'
