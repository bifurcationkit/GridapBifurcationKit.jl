# [Simple Hopf point](@id simple-hopf)

At a Hopf branch point ``(x_0,p_0)`` for the problem ``F(x,p)=0``, the spectrum of the linear operator ``dF(x_0,p_0)`` contains two purely imaginary eigenvalues ``\pm i\omega,\ \omega > 0`` which are simple. At such a point, we can compute the **normal form** to transform the problem

```math
\dot x = \mathbf{F}(x;p)
```

in large dimensions to a **complex** polynomial vector field (``\delta p\equiv p-p_0``):

```math
\dot z = z\left(a \cdot\delta p + i\omega + l_1|z|^2\right)\quad\text{(E)}
```

whose solutions give access to the solutions of the Cauchy problem in a neighborhood of ``(x,p)``.

!!! tip "Coefficient $l_1$"
    The coefficient ``l_1`` above is called the **Lyapunov** coefficient.

!!! note "Differential-algebraic problems"
    For a problem with a mass matrix (a `GridapBifProblem`), the Hopf point is characterized
    by a pair of purely imaginary eigenvalues of the **pencil** ``(dF, M)``, see
    [Mass matrix](@ref mass-matrix).

## Normal form computation

The normal form (E) is automatically computed as follows

```julia
get_normal_form(br::ContResult, ind_bif::Int ;
	verbose = false, ζs = nothing, lens = br.param_lens)
```

where `br` is a branch computed after a call to [`continuation`](@ref) with detection of
bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch
`br`. The above call returns a point with information needed to compute the bifurcated
branch. For more information about the optional parameters, we refer to
[`get_normal_form`](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/).

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

## See also

- [Branch switching](@ref Branch-switching-page)
- [Fold / Hopf Continuation](@ref)

## References

[^Haragus]: > Haragus, Mariana, and Gérard Iooss. Local Bifurcations, Center Manifolds, and Normal Forms in Infinite-Dimensional Dynamical Systems. London: Springer London, 2011. https://doi.org/10.1007/978-0-85729-112-7.
