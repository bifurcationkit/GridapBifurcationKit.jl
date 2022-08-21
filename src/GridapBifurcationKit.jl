module GridapBifurcationKit
	using Gridap, BifurcationKit, Setfield
	const BK = BifurcationKit

	include("Wrap.jl")
	include("LinearSolvers.jl")
	include("Bifurcation.jl")

	export GridapBifProblem, reactionDiffusion
	export isSymmetric, residual, jacobian, d2F, d3F
	export newton, continuation, computeNormalForm, bifurcationdiagram

end # module
