# GridapBifurcationKit.jl

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/mittelmannGridap/#d-Bratu%E2%80%93Gelfand-problem-with-[Gridap.jl](https://github.com/gridap/Gridap.jl)-(Intermediate)) |  |


This Julia package aims at performing **automatic bifurcation analysis** of PDE solved using the Finite Elements Method (FEM) with the Julia package [Gridap.jl](https://github.com/gridap/Gridap.jl).

> I would like to thank **Santiago Badia** and **Francesc Verdugo** for their help in developing this package.



**If you use this package for your work, please cite it!! Open source development strongly depends on this. It is referenced as follows:**

```
@misc{veltz:hal-02902346,
  TITLE = {{BifurcationKit.jl}},
  AUTHOR = {Veltz, Romain},
  URL = {https://hal.archives-ouvertes.fr/hal-02902346},
  INSTITUTION = {{Inria Sophia-Antipolis}},
  YEAR = {2020},
  MONTH = Jul,
  KEYWORDS = {pseudo-arclength-continuation ; periodic-orbits ; floquet ; gpu ; bifurcation-diagram ; deflation ; newton-krylov},
  PDF = {https://hal.archives-ouvertes.fr/hal-02902346/file/354c9fb0d148262405609eed2cb7927818706f1f.tar.gz},
  HAL_ID = {hal-02902346},
  HAL_VERSION = {v1},
}

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

## Installation 

To install this package, run the command

```julia
add https://github.com/bifurcationkit/GridapBifurcationKit.jl
```


## Main features

Most [features](https://github.com/rveltz/BifurcationKit.jl#main-features) of [BifurcationKit](https://github.com/rveltz/BifurcationKit.jl) are ported, please see the `examples` folder or the [tutorials](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/tutorials/mittelmannGridap/#d-Bratu–Gelfand-problem-with-[Gridap.jl](https://github.com/gridap/Gridap.jl)-(Intermediate)) for example of use.
