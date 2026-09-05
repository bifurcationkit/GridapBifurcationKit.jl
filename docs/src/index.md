# GridapBifurcationKit.jl

This Julia package aims at performing **automatic bifurcation analysis of partial differential equations** (PDEs) discretized with the finite element toolbox [Gridap.jl](https://github.com/gridap/Gridap.jl).

It builds upon [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl) to perform continuation and numerical bifurcation analysis. The equations are written in a weak form very close to the mathematical formulation, and all the (large scale) Newton-Krylov machinery of `BifurcationKit.jl` is available.

## 📦 Installation

Assuming that you already have Julia correctly installed, it suffices to import `GridapBifurcationKit.jl` in the standard way:

```julia
import Pkg
Pkg.add("https://github.com/bifurcationkit/GridapBifurcationKit.jl")
```

## Main features

- Newton-Krylov solver with generic linear / eigen *preconditioned* solver, and arc-length continuation.
- Continuation methods: PALC, Moore-Penrose, etc. See [methods](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/IntroContinuation/).
- Monitoring user functions along curves computed by continuation, see [events](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/EventCallback/).
- Bifurcation points are located using a bisection algorithm.
- Detection of Branch, Fold, Hopf bifurcation points of stationary solutions and computation of their normal form.
- Automatic branch switching at simple Hopf points to periodic orbits.
- **Automatic bifurcation diagram computation of equilibria.**
- Fold / Hopf continuation.
- Support for **differential-algebraic problems** through a mass matrix: PDEs of the form

```math
M\frac{\mathrm{d}z}{\mathrm{d}t} = F(z, p),
```

with a possibly singular ``M`` (for example incompressible flows, where the pressure has no time derivative). This is essential to correctly compute the stability of the solution and to detect Hopf bifurcations.

Custom state means, we can use something else than `AbstractArray`, for example your own `struct`.

|Features| DAE (mass)| Matrix Free|Custom state| Tutorials |
|---|---|---|---|---|
| (Deflated) Krylov-Newton| Yes | Yes | Yes| |
| Continuation PALC (Natural, Secant, Tangent, Polynomial) | Yes | | | |
| Bifurcation / Fold / Hopf point detection | Yes | Y | | |
| Fold / Hopf point continuation | Yes | Y | | |
| Branch switching at Branch / Hopf points | Yes | Y | `AbstractArray` | |
| Automatic bifurcation diagram computation of equilibria | Yes | Y | `AbstractArray` | |
| Periodic orbit (collocation / trapeze) Newton / continuation | Yes | | `AbstractVector` | |
| Codim 2 bifurcation detection (BT, GH, cusp, ZH, HH) | Yes | Y | | |

## 🧑‍💻 Other softwares

There are several good softwares already available for PDE bifurcation analysis. One can mention e.g. [`pde2path`](https://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/). The present package is, to our knowledge, the only one written in Julia and leveraging automatic finite element discretization with `Gridap.jl`.

## 📚 Citing this work

If you use this package for your work, we ask that you **cite** the following paper! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).

You need to cite this entry **as well**:

```
@article{Badia2020,
  doi = {10.21105/joss.02520},
  url = {https://doi.org/10.21105/joss.02520},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2520},
  author = {Santiago Badia and Francesc Verdugo},
  title = {Gridap: An extensible Finite Element toolbox in Julia},
  journal = {Journal of Open Source Software}
}
```

## Citations

Papers citing this work are collected on [Zotero](https://www.zotero.org/groups/6097154/citations_of_bifurcationkit/library).

These citations are aggregated from [Google Scholar (search)](https://scholar.google.com/scholar?q=bifurcationkit&hl=en&as_sdt=0,5) and [Google Scholar (citations)](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=159498619004863176,12573642401780006854,8662907770106865595).
