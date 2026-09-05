# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Parameters

```@docs
BifurcationKit.NewtonPar
```

```@docs
BifurcationKit.ContinuationPar
```

## Problems

```@docs
GridapBifProblem
```

```@docs
GridapBifurcationKit.get_mass_matrix
```

## Eigen solvers

The generic eigensolvers are provided by `BifurcationKit.jl` (`DefaultEig`, `EigArpack`,
`EigArnoldiMethod`, `EigKrylovKit`) and the DAE wrappers by `EigenDAE` and
`EigenMassMatrix`. See [Eigen Solvers](@ref eigensolver-page).

```@docs
BifurcationKit.EigenMassMatrix
```

## Branch switching (branch point)

Automatic branch switching at a branch point is performed with

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

where `br` is a branch computed with detection of bifurcation points enabled. See
[Branch switching](@ref Branch-switching-page) for more information and the precise method
definition in `BifurcationKit.jl`.

## Branch switching (Hopf point)

Automatic branch switching at a Hopf point towards periodic orbits is performed with

```julia
continuation(br, ind_bif::Int, _contParams::ContinuationPar,
	prob::AbstractPeriodicOrbitProblem ; δp = nothing, ampfactor = 1, kwargs...)
```

See [Branch switching](@ref Branch-switching-page) for more information and the precise
method definition in `BifurcationKit.jl`.

## Normal form

The normal forms of simple branch points and Hopf points are computed with the
`BifurcationKit` function `get_normal_form`, see [Simple bifurcation branch point](@ref simple-bp) and
[Simple Hopf point](@ref).
