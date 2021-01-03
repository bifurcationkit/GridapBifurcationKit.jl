using Revise
using Plots
using Gridap,Setfield
using Gridap.FESpaces
using GridapBifurcationKit, BifurcationKit
plotgridap!(x; k...) = (n=Int(sqrt(length(x)));heatmap!(reshape(x,n,n); color=:viridis, k...))
plotgridap(x; k...) =( plot();plotgridap!(x; k...))
norminf = x -> norm(x, Inf)
#############################################
# discretisation
n = 40
domain = (0,1,0,1)
cells = (n,n)
model = CartesianDiscreteModel(domain,cells)

# function spaces
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,conformity=:H1,)#dirichlet_tags="boundary")
U = TrialFESpace(V)

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)


NL(u) = exp(u)
res(u,p,v) = ∫( -∇(v)⋅∇(u) -  v ⋅ (u - p.λ ⋅ (NL ∘ u)) * 10 )*dΩ
jac(u,p,du,v) = ∫( -∇(v)⋅∇(du) - v ⋅ du ⋅ (1 - p.λ *( NL ∘ u)) * 10 )*dΩ
d2res(u,p,du1,du2,v) = ∫( v ⋅ du1 ⋅ du2 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ
d3res(u,p,du1,du2,du3,v) = ∫( v ⋅ du1 ⋅ du2 ⋅ du3 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ

uh = zero(U)
par_bratu=(λ = 0.01,)
# prob = GridapProblem(res, jac, d2res, d3res, V, U)
prob = GridapProblem(res, V, U)

optn = NewtonPar(eigsolver = EigArpack())#EigKrylovKit(dim = 100))
sol, = newton(prob, uh, par_bratu, NewtonPar(optn; verbose = true))

# weight for normbratu
const w = cumsum(ones(length(uh.free_values))) / length(uh.free_values)
w .= (1 .+ LinRange(-1,1,n+1)) * (LinRange(-1,1,n+1))' |> vec
w .-= minimum(w)
normbratu = x -> norm(x .* w) / sqrt(length(x))

opts = ContinuationPar(pMax = 40., pMin = 0.01, ds = 0.01, maxSteps = 1000, detectBifurcation = 3, newtonOptions = optn, nev = 20, precisionStability = 1e-6, nInversion = 6, dsminBisection = 1e-17, maxBisectionSteps = 25, tolBisectionEigenvalue = 1e-19)
br, = continuation(prob, uh, par_bratu, (@lens _.λ), opts;
	plot = true,
	verbosity = 0,
	plotSolution = (x,p;k...) -> plotgridap!(x;  k...),
	# printSolution = (x, p) -> normbratu(x),
	# finaliseSolution = finSol,
	)

plot(br)

nf = computeNormalForm(prob, br, 3; verbose = true, nev = 50)
####################################################################################################
br1, = continuation(prob, br, 3,
		setproperties(opts;ds = 0.005, dsmax = 0.05, maxSteps = 140, detectBifurcation = 3);
		verbosity = 0, plot = true, nev = 10,
		# printSolution = (x, p) -> normbratu(x),
		# finaliseSolution = finSol,
		tangentAlgo = BorderedPred(),
		usedeflation = true,
		# callbackN = cb,
		plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
		)

plot(br,br1)

br2, = continuation(prob, br1, 3,
		setproperties(opts;ds = 0.005, dsmax = 0.05, maxSteps = 140, detectBifurcation = 0);
		verbosity = 0, plot = true, nev = 10,
		# printSolution = (x, p) -> normbratu(x),
		# finaliseSolution = finSol,
		tangentAlgo = BorderedPred(),
		usedeflation = true,
		# callbackN = cb,
		plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
		)

plot(br,br1,br2, legend=false)

br3, = continuation(prob, br, 2,
		setproperties(opts;ds = 0.005, dsmax = 0.05, maxSteps = 140, detectBifurcation = 3);
		verbosity = 0, plot = true, nev = 10,
		# printSolution = (x, p) -> normbratu(x),
		# finaliseSolution = finSol,
		# tangentAlgo = BorderedPred(),
		usedeflation = true,
		# callbackN = cb,
		plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
		)

plot(br,br1,br2,br3..., legend=false)
plot(br3)
####################################################################################################
function optionsCont(x,p,l; opt = opts)
	if l <= 1
		return opt
	elseif l==2
		@info "l2"
		return setproperties(opt ;detectBifurcation = 3,ds = 0.005, dsmax = 0.025, a = 0.75)
	else
		@info "l>2"
		return setproperties(opt ;detectBifurcation = 0,ds = 0.00051, dsmax = 0.05, dsmin =  0.0005)
	end
end

function cb(x,f,J,res,it,itl,optN; kwargs...)
	fromnewton = get(kwargs, :fromNewton, false)
	if ~fromnewton
		_x = get(kwargs, :z0, nothing)
		cdt = (norm(_x.u - x) < 1020.5 && abs(_x.p - kwargs[:p]) < 0.05)
		~cdt && @warn "Jump callback"
		return cdt
	end
	true
end

function finSol(z, tau, step, br)
	if ~isnothing(br.bifpoint)
		if br.bifpoint[end].step == step
			BifurcationKit._show(stdout, br.bifpoint[end], step)
		end
	end
	return true
end

diagram = GridapBifurcationKit.bifurcationdiagram(prob,
		uh, par_bratu, (@lens _.λ), 2, optionsCont;
		# δp = 0.001,
		verbosity = 0, plot = true,
		printSolution = (x, p) -> normbratu(x),
		callbackN = cb,
		# tangentAlgo = BorderedPred(),
		usedeflation = true,
		finaliseSolution = finSol,
		plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
		# normC = norminf
		)

code = (2,)
	plot(diagram; code = code, level = (0, 5), plotfold = false, putbifptlegend=false, markersize=2)
	plot!(br)
	# xlims!(0.01, 0.4)
	title!("#branches = $(size(getBranch(diagram, code)))")
	# xlims!(0.01, 0.065, ylims=(2.5,6.5))
plot!(br2)

diagram1 = GridapBifurcationKit.bifurcationdiagram(prob,
	# this improves the first branch on the violet curve. Note that
	# for symmetry reasons, the first bifurcation point
	# has 8 branches
	br, 3, optionsCont;
	verbosity = 0, plot = true,
	printSolution = (x, p) -> normbratu(x),
	callbackN = cb,
	finaliseSolution = finSol,
	usedeflation = true,
	plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
	normC = norminf,
	)

GridapBifurcationKit.bifurcationdiagram!(prob,
	getBranch(diagram, (2,)), (current = 3, maxlevel = 6), optionsCont;
	verbosity = 0, plot = true,
	printSolution = (x, p) -> normbratu(x),
	callbackN = cb,
	finaliseSolution = finSol,
	usedeflation = true,
	plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
	normC = norminf,
	)

diagram1 = GridapBifurcationKit.bifurcationdiagram(prob,
	# this improves the first branch on the violet curve. Note that
	# for symmetry reasons, the first bifurcation point
	# has 8 branches
	getBranch(diagram, (3,)).γ, 2, optionsCont;
	verbosity = 0, plot = true,
	# printSolution = (x, p) -> normbratu(x),
	callbackN = cb,
	finaliseSolution = finSol,
	usedeflation = true,
	plotSolution = (x, p; k...) -> plotgridap!(x;  k...),
	normC = norminf,
	)
