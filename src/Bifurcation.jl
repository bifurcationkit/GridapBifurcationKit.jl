import BifurcationKit: newton, continuation, computeNormalForm

# regular newton method
function newton(prob::GridapProblem, x0, par, options::NewtonPar; kwargs...)
	BK.newton(
		(u, p) -> prob(Val(:Res), u, p),
		(u, p) -> prob(Val(:Jac), u, p),
		get_free_dof_values(x0), par, options; kwargs...)
end

# simple continuation
function continuation(prob::GridapProblem, x0, par, lens::Lens, contParams::ContinuationPar; linearAlgo = nothing, kwargs...)
	BK.continuation(
	(u, p) -> prob(Val(:Res), u, p),
	(u, p) -> prob(Val(:Jac), u, p),
	get_free_dof_values(x0), par, lens, contParams; kwargs...)
end

# normal form computation
function computeNormalForm(prob::GridapProblem, br::BK.AbstractBranchResult, ind_bif::Int; kwargs...)
	BK.computeNormalForm(
		(u, p) -> prob(Val(:Res), u, p),
		(u, p) -> prob(Val(:Jac), u, p),
		(u, p, dx1, dx2) -> prob(u, p, dx1, dx2),
		(u, p, dx1, dx2, dx3) -> prob(u, p, dx1, dx2, dx3),
		br, ind_bif; kwargs...
	)
end

# automatic branch switching
function continuation(prob::GridapProblem, br::BK.AbstractBranchResult, ind_bif::Int, contParams::ContinuationPar; kwargs...)
	BK.continuation(
		(u, p) -> prob(Val(:Res), u, p),
		(u, p) -> prob(Val(:Jac), u, p),
		(u, p, dx1, dx2) -> prob(u, p, dx1, dx2),
		(u, p, dx1, dx2, dx3) -> prob(u, p, dx1, dx2, dx3),
		br, ind_bif, contParams; kwargs...
	)
end

# automatic bifurcation diagram
function bifurcationdiagram(prob::GridapProblem, x0, par0, lens::Lens, level::Int, options; usedeflation = false, kwargs...)
	BK.bifurcationdiagram(
		(u, p) -> prob(Val(:Res), u, p),
		(u, p) -> prob(Val(:Jac), u, p),
		(u, p, dx1, dx2) -> prob(u, p, dx1, dx2),
		(u, p, dx1, dx2, dx3) -> prob(u, p, dx1, dx2, dx3),
		get_free_dof_values(x0), par0, lens, level, options; usedeflation = usedeflation, kwargs...)
end

function bifurcationdiagram(prob::GridapProblem, br::BK.AbstractBranchResult, level::Int, options; usedeflation = false, kwargs...)
	BK.bifurcationdiagram(
		(u, p) -> prob(Val(:Res), u, p),
		(u, p) -> prob(Val(:Jac), u, p),
		(u, p, dx1, dx2) -> prob(u, p, dx1, dx2),
		(u, p, dx1, dx2, dx3) -> prob(u, p, dx1, dx2, dx3),
		br, level, options; usedeflation = usedeflation, kwargs...)
end

function bifurcationdiagram!(prob::GridapProblem, node::BK.BifDiagNode, level::NamedTuple{(:current, :maxlevel), Tuple{Int64,Int64}}, options; code = "0", usedeflation = false, kwargs...)
	BK.bifurcationdiagram!(
		(u, p) -> prob(Val(:Res), u, p),
		(u, p) -> prob(Val(:Jac), u, p),
		(u, p, dx1, dx2) -> prob(u, p, dx1, dx2),
		(u, p, dx1, dx2, dx3) -> prob(u, p, dx1, dx2, dx3),
		node, level, options; code = code, usedeflation = usedeflation, kwargs...)
end
