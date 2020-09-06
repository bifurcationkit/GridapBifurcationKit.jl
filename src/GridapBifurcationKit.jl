module GridapBifurcationKit
	using Gridap, BifurcationKit, Setfield
	const BK = BifurcationKit

	include("Wrap.jl")
	include("LinearSolvers.jl")
	include("Bifurcation.jl")

	export GridapProblem
	export newton, continuation, computeNormalForm, bifurcationdiagram

end # module
