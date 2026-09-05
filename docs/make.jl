using Pkg
cd(@__DIR__)
pkg" activate ."

# pkg"dev Gridap"
# pkg"dev /Users/rveltz/work/prog_gd/julia/dev/dev1/bkorg/GridapBifurcationKit"

using Documenter, GridapBifurcationKit, BifurcationKit
using DocumenterCodeBlocks
ENV["GKSwstype"] = "100"

# to display progress
ENV["JULIA_DEBUG"] = Documenter

format = Documenter.HTML(;
		collapselevel = 1,
		size_threshold_warn = 300 * 2^10, # raise slightly from 100 to 200 KiB
		size_threshold = 400 * 2^10,      # raise slightly 200 to to 300 KiB
		assets=[
			asset("https://bifurcationkit.github.io/assets/js/documentation.js"),
			asset("https://bifurcationkit.github.io/assets/css/documentation.css"),
				],
		)# assets = ["assets/indigo.css"]),

makedocs(
	modules = [GridapBifurcationKit, BifurcationKit],
	doctest = false,
	pagesonly = false, # this is on Documenter#master, do not compile what is not in pages =
	draft = false,
	warnonly = true,
	sitename = "Bifurcation of PDEs based on Gridap in Julia",
	format = format,
	authors = "Romain Veltz",
	plugins = [CodeBlocks()],
	pages = Any[
		"🏠 Home" => "index.md",
		"📎 Tutorials" => "tutorials/tutorials.md",
		"🧩 Problems" => [
			"Bifurcation Problem" => "problems.md",
		],
		"🧰 Functionalities" => [
			"Bifurcations" => [
				"Bifurcation detection (codim 1)" => "detectionBifurcation.md",
				"Fold / Hopf Continuation (codim 2)" => "codim2Continuation.md",
			],
			"Normal form" => [
				"Simple branch point" => "simplebp.md",
				"Simple Hopf point" => "simplehopf.md",
			],
			"Branch switching" => "branchswitching.md",
		],
		"⚙️ Options" => [
			"Eigen Solvers" => "eigensolver.md",
		],
		"❓ Frequently Asked Questions" => "faq.md",
		"📚 Library" => "library.md",
	]
	)

deploydocs(
	repo = "github.com/bifurcationkit/GridapBifurcationKit.jl.git",
	push_preview = true, 
	target = "build", 
	devbranch = "master"
	)
