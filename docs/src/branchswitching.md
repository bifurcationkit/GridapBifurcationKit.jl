# [Branch switching](@id Branch-switching-page)

The precise definition of the methods is given in
[Branch switching (branch point)](@ref) and
[Branch switching (Hopf point)](@ref).

```@contents
Pages = ["branchswitching.md"]
Depth = 3
```

## Branch switching from simple branch point to equilibria

You can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

where `br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled. This call computes the branch bifurcating from the `ind_bif`-th bifurcation point in `br`. An example of use is provided in the [1d Bratu model](@ref bratu) tutorial.

## Branch switching from non-simple branch point to equilibria

We provide an automatic branch switching method in this case. The method is to first compute the reduced equation (see [Simple bifurcation branch point](@ref simple-bp)) and use it to compute the nearby solutions. These solutions are seeded as initial guess for [`continuation`](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

This requires the second and third derivatives of the residual (`d2res`, `d3res`) to be provided when building the [`GridapBifProblem`](@ref).

## Branch switching from Hopf point to periodic orbits

In order to compute the bifurcated branch of periodic solutions at a Hopf bifurcation point, you need to choose a method to compute periodic orbits among:

- periodic orbits based on orthogonal collocation (`PeriodicOrbitOCollProblem`),
- periodic orbits based on the trapezoidal rule (`Trapeze`).

Once you have decided which method to use, you use the following call:

```julia
continuation(br::ContResult, ind_HOPF::Int, _contParams::ContinuationPar,
	prob::AbstractPeriodicOrbitProblem ;
	δp = nothing, ampfactor = 1, kwargs...)
```

We refer to [`continuation`](@ref) for more information about the arguments.

!!! note "Mass matrix and periodic orbits"
    For a problem with a mass matrix (see [Mass matrix](@ref mass-matrix)), one builds the
    periodic orbit problem with `massmatrix = M`, for example

    ```julia
    probPO = Trapeze(; M, N = length(u0), massmatrix = GridapBifurcationKit.get_mass_matrix(prob))
    ```

    This accounts for the differential-algebraic structure ``M\dot z = F(z, p)``.
