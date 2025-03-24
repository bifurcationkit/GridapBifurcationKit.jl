# GridapBifurcationKit.jl

This Julia package aims at performing **bifurcation analysis** of Homoclinic / Heteroclinic orbits of Cauchy problems.

It builds upon [BifurcationKit.jl]() with version > 0.2 to perform continuation and numerical bifurcation analysis.

## 📦 Installation

Assuming that you already have Julia correctly installed, it suffices to import  `GridapBifurcationKit.jl` in the standard way:

`import Pkg; Pkg.add("https://github.com/bifurcationkit/GridapBifurcationKit.jl")`

## 📚 Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).
```

You need to cite this entry **as well**

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