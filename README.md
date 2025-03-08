# GridapBifurcationKit.jl

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/mittelmannGridap/#d-Bratu%E2%80%93Gelfand-problem-with-[Gridap.jl](https://github.com/gridap/Gridap.jl)-(Intermediate)) |  |


This Julia package aims at performing **automatic bifurcation analysis** of PDE solved using the Finite Elements Method (FEM) with the Julia package [Gridap.jl](https://github.com/gridap/Gridap.jl).

> I would like to thank **Santiago Badia** and **Francesc Verdugo** for their help in developing this package.

## 📦 Installation 

To install this package, run the command

```julia
add https://github.com/bifurcationkit/GridapBifurcationKit.jl
```

## 📚 Support and citation
If you use `BifurcationKit.jl` in your work, we ask that you cite the following paper on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib). **Open source development strongly depends on this.**

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


## Main features

Most [features](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/capabilities/) of [BifurcationKit](https://github.com/rveltz/BifurcationKit.jl) are ported, please see the `examples` folder or the [tutorials](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/tutorials/mittelmannGridap/#d-Bratu–Gelfand-problem-with-[Gridap.jl](https://github.com/gridap/Gridap.jl)-(Intermediate)) for example of use.
