module GridapBifurcationKit
	using Gridap, BifurcationKit, Setfield
	using FiniteDifferences
	const BK = BifurcationKit

	include("Wrap.jl")
	include("LinearSolvers.jl")
	include("Bifurcation.jl")
	include("reactionDiffusion.jl")

	export GridapBifProblem, reactionDiffusion
	export isSymmetric, residual, jacobian, d2F, d3F
	export newton, continuation, computeNormalForm, bifurcationdiagram

end # module
