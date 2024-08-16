using Pkg
cd(@__DIR__)
pkg" activate ."

pkg"dev Gridap"
pkg"dev /Users/rveltz/work/prog_gd/julia/dev/dev1/bkorg/GridapBifurcationKit"

using Documenter, GridapBifurcationKit, Setfield, BifurcationKit
ENV["GKSwstype"] = "100"

# to display progress
ENV["JULIA_DEBUG"] = Documenter

makedocs(	
	modules = [BifurcationKit],
	doctest = false,
	pagesonly = false, # this is on Documenter#master, do not compile what is not in pages =
	draft = false,
	warnonly = true,
	sitename = "Bifurcation of PDEs based on Gridap in Julia",
	format = Documenter.HTML(collapselevel = 1,),# assets = ["assets/indigo.css"]),
	authors = "Romain Veltz",
	pages = Any[
		"Home" => "index.md",
		"Tutorials" => "tutorials/tutorials.md",
		"Problems" => [
			"bifurcation problem" => "problems.md",
		],
		"Functionalities" => [
			"Bifurcations" => [
				"Bifurcation detection (codim 1)" => "detectionBifurcation.md",
				"Branch switching" => "branchswitching.md",
							  ],		
		],
		"Frequently Asked Questions" => "faq.md",
		"Library" => "library.md"
	]
	)

deploydocs(
	repo = "github.com/bifurcationkit/GridapBifurcationKit.jl.git",
	devbranch = "main"
)
