import BifurcationKit: _getvectortype

# Mass operator flavours. `mass_type` is inferred from the `mass` keyword of
# `GridapBifProblem` and drives `get_mass_matrix` dispatch:
#   * `MassDefaut()`: no mass given, use the L² mass ∫(u⋅v)dΩ;
#   * `ConstMass()` : a constant (state-independent) integrand `mass(u, v)`;
#   * `SDMass()`    : a state-dependent operator `mass(u, p, du, v)`, i.e. M(x, p).
abstract type AbstractMassType end
struct MassDefaut <: AbstractMassType end
struct ConstMass <: AbstractMassType end
struct SDMass <: AbstractMassType end

# infer the mass type from a user-provided `mass` keyword: a 4-argument method
# `mass(u, p, du, v)` is reported by `methods` as `nargs == 5` (callable + 4).
_mass_type(::Nothing) = MassDefaut()
_mass_type(mass) = any(m -> m.nargs == 5, methods(mass)) ? SDMass() : ConstMass()

_is_constant_mass(::MassDefaut) = true
_is_constant_mass(::ConstMass) = true
_is_constant_mass(::SDMass) = false

# Analytic parameter-derivative forms. A `R01`/`R02` residual-like form has
# signature `(u, p, v)` (callable + 3 args = `nargs == 4`); a `R11` jacobian-like
# form has signature `(u, p, du, v)` (callable + 4 args = `nargs == 5`). Anything
# else (BifurcationKit sentinels `FiniteDifferences()`/`AutoDiff()`, `nothing`, or
# an already assembled `(x, p) -> Vector` closure) is forwarded untouched.
_is_residual_form(f) = f isa Function && any(m -> m.nargs == 4, methods(f))
_is_jac_form(f)      = f isa Function && any(m -> m.nargs == 5, methods(f))

"""
    GridapProblem(res, jac, d2res, d3res, V, U, ls, dΩ, mass, mass_type)

Low-level, **parameter-agnostic** description of a PDE problem discretized with
`Gridap`. It is the object stored in the `probFE` field of a
[`GridapBifProblem`](@ref) and is not meant to be built directly by the user:
use [`GridapBifProblem`](@ref) instead.

Contrary to [`GridapBifProblem`](@ref), the parameters are **not** stored here:
the forms take an explicit parameter argument `p`. This keeps the object
independent of the parameter value, so that it can be reused (and cheaply
re-wrapped in an `FEOperator` by `op_from_param`) along a continuation.

# Fields
- `res(u, p, v)`: residual of the (semi-)discretized problem. It returns a
  `Gridap` `DomainContribution` and is the (weak) form whose zero is sought.
- `jac(u, p, du, v)`: jacobian of `res` with respect to the state, in the
  direction `du`. If `nothing`, the jacobian is obtained by finite differences
  when the corresponding `FEOperator` is built.
- `d2res(u, p, du1, du2, v)`: second derivative (Hessian) of `res`. May be
  `nothing`, in which case it is approximated by finite differences.
- `d3res(u, p, du1, du2, du3, v)`: third derivative (Tressian) of `res`. May be
  `nothing`, in which case it is approximated by finite differences.
- `V`: the `TestFESpace`.
- `U`: the `TrialFESpace` (it carries the Dirichlet data).
- `ls`: reserved linear solver, or `nothing` (currently unused).
- `dΩ`: the `Measure` used to assemble the residual; it is stored so that the
  mass matrix can be assembled later (see [`get_mass_matrix`](@ref)).
- `mass`: mass operator, or `nothing` for the default L² mass. See
  [`get_mass_matrix`](@ref) for the accepted signatures.
- `mass_type`: kind of mass operator, one of `MassDefaut()` (default L² mass),
  `ConstMass()` (state-independent bilinear integrand) or `SDMass()`
  (state-dependent `mass(u, p, du, v)`). It is inferred from `mass` by
  `_mass_type`.

# Callable / extended methods
- `(gp::GridapProblem)(u, p, du1, du2)` returns the vector associated to
  `d2res(u, p, du1, du2, v)` (used by `BifurcationKit.d2F`).
- `(gp::GridapProblem)(u, p, du1, du2, du3)` returns the vector associated to
  `d3res(u, p, du1, du2, du3, v)` (used by `BifurcationKit.d3F`).
- `residual` and `jacobian` evaluate `res` and `jac` at a given state and
  parameter value.
- [`get_mass_matrix`](@ref) assembles the mass matrix.

# See also
- [`GridapBifProblem`](@ref), [`get_mass_matrix`](@ref)
"""
struct GridapProblem{Tres, Tjac, Td2res, Td3res, TV, TU, Tls, TOm, Tm, Tmt}
    res::Tres        # res(u, p, v),                 residual
    jac::Tjac        # jac(u, p, du, v),             jacobian
    d2res::Td2res    # d2res(u, p, du1, du2, v)
    d3res::Td3res    # d3res(u, p, du1, du2, du3, v)
    V::TV
    U::TU
    ls::Tls
    dΩ::TOm
    mass::Tm
    mass_type::Tmt
end

_deriv_name(f) = f === nothing ? "missing" : "provided"

function Base.show(io::IO, gp::GridapProblem; prefix = "")
    print(io, prefix, "Gridap FE problem\n")
    print(io, prefix, "├─ Trial/Test:  ", nameof(typeof(gp.U)), " / ", nameof(typeof(gp.V)), "\n")
    print(io, prefix, "├─ Free dofs:   ", Gridap.FESpaces.num_free_dofs(gp.U), "\n")
    print(io, prefix, "├─ Mass:        ", nameof(typeof(gp.mass_type)), "\n")
    print(io, prefix, "└─ Derivatives: jac=", _deriv_name(gp.jac),
          ", d2res=", _deriv_name(gp.d2res), ", d3res=", _deriv_name(gp.d3res))
end

# rebuild a gridap operator for each parameter value
function op_from_param(gp::GridapProblem{Tres, Tjac}, p) where {Tres, Tjac}
    res(u, v) = gp.res(u, p, v)
    jac(u, du, v) = gp.jac(u, p, du, v)
    return FEOperator(res, jac, gp.U, gp.V)
end

function op_from_param(gp::GridapProblem{Tres, Nothing}, p) where {Tres}
    res(u, v) = gp.res(u, p, v)
    return FEOperator(res, gp.U, gp.V)
end

# residual
function residual(gp::GridapProblem, u::AbstractArray{ <: Real}, p)
    op = op_from_param(gp, p)
    algop = Gridap.FESpaces.get_algebraic_operator(op)
    return Gridap.FESpaces.residual(algop, u)
end

# (sparse) jacobian matrix
function jacobian(gp::GridapProblem, u, p)
    op = op_from_param(gp, p)
    algop = Gridap.FESpaces.get_algebraic_operator(op)
    return Gridap.FESpaces.jacobian(algop, u)
end

# assemble an analytic parameter-derivative residual form `form(u, p, v)` into a
# free dof vector (same convention as `residual(gp, u, p)`)
function _residual_from_form(gp::GridapProblem, form, x, p)
    op   = FEOperator((u, v) -> form(u, p, v), gp.U, gp.V)
    alop = Gridap.FESpaces.get_algebraic_operator(op)
    return Gridap.FESpaces.residual(alop, x)
end

# assemble an analytic parameter-derivative jacobian form `form(u, p, du, v)` and
# apply it to `dx` (used for `R11`)
function _apply_from_biform(gp::GridapProblem, form, x, p, dx)
    uh = FEFunction(gp.U, x)
    A  = Gridap.FESpaces.assemble_matrix((du, v) -> form(uh, p, du, v), gp.U, gp.V)
    return A * dx
end

# second derivative
function (gp::GridapProblem)(u, p, du1, du2)
    du1h = FEFunction(gp.U, du1)
    du2h = FEFunction(gp.U, du2)
    a(u, v) = gp.d2res(u, p, du1h, du2h, v)
    feop = FEOperator(a, gp.U, gp.V)
    alop = Gridap.FESpaces.get_algebraic_operator(feop)
    Gridap.FESpaces.residual(alop, u)
end

# third derivative
function (gp::GridapProblem)(u, p, du1, du2, du3)
    du1h = FEFunction(gp.U, du1)
    du2h = FEFunction(gp.U, du2)
    du3h = FEFunction(gp.U, du3)
    a(u, v) = gp.d3res(u, p, du1h, du2h, du3h, v)
    feop = FEOperator(a, gp.U, gp.V)
    alop = Gridap.FESpaces.get_algebraic_operator(feop)
    Gridap.FESpaces.residual(alop, u)
end

# second derivative
function (gp::GridapProblem{Tres, Tjac, Nothing})(u, p, du1, du2) where {Tres, Tjac}
    jvp(central_fdm(3, 1), z -> jacobian(gp, z, p) * du1, (u, du2))
end

# third derivative
function (gp::GridapProblem{Tres, Tjac, Td2res, Nothing})(u, p, du1, du2, du3) where {Tres, Tjac, Td2res}
    jvp(central_fdm(3, 1), z -> gp(z, p, du1, du2), (u, du3))
end

mass_default(u,v,dΩ) = ∫(u⋅v) * dΩ

# shared assembly of a (bi)linear mass integrand over the free dofs
function _assemble_mass_matrix(gp::GridapProblem, integrand)
    (;U, V) = gp
    v = get_fe_basis(V)
    u = get_trial_fe_basis(U)
    assemblytuple = Gridap.FESpaces.collect_cell_matrix(U,V,integrand(u,v))
    cell_matrix_MM   = collect(assemblytuple[1][1]) # This result is no longer a LazyArray
    newassemblytuple = ([cell_matrix_MM], assemblytuple[2], assemblytuple[3])
    a = SparseMatrixAssembler(U, V)
    return assemble_matrix(a, newassemblytuple)
end

"""
    get_mass_matrix(gp::GridapBifProblem[, x, p])
    get_mass_matrix(gp::GridapProblem[, x, p], dΩ = gp.dΩ)

Assemble the (sparse) mass matrix associated to the problem, on the free dofs of the
trial/test spaces (Dirichlet dofs are eliminated, consistently with the jacobian).

The operator is selected by `gp.mass_type`, itself inferred from the `mass`
keyword of [`GridapBifProblem`](@ref):

* `MassDefaut()`: no mass given, use the L² mass `∫(u⋅v)*dΩ`;
* `ConstMass()` : the constant (state-independent) integrand `mass(u, v)`;
* `SDMass()`    : the state-dependent operator `mass(u, p, du, v)`, reconstructed
  from the free dof vector `x` (i.e. `M(x, p)`); the `(gp, x, p)` methods must
  be used.

For an incompressible flow, a typical constant choice is the *velocity only*
mass `(u,p),(v,q) -> ∫(v⊙u)*dΩ`, which yields a singular mass matrix with a zero
pressure block, as required for the stability of the differential-algebraic
system ``M\\dot z = F(z, p)``.
"""
function get_mass_matrix(gp::GridapProblem, dΩ = gp.dΩ)
    return _get_mass_matrix(gp.mass_type, gp, dΩ)
end

_get_mass_matrix(::MassDefaut, gp::GridapProblem, dΩ) =
    _assemble_mass_matrix(gp, (u, v) -> mass_default(u, v, dΩ))

_get_mass_matrix(::ConstMass, gp::GridapProblem, dΩ) =
    _assemble_mass_matrix(gp, gp.mass)

function _get_mass_matrix(::SDMass, ::GridapProblem, dΩ)
    throw(ArgumentError("The mass operator passed to `GridapBifProblem` is state-dependent; call `get_mass_matrix(gp, x, p)` with the state `x` and parameters `p`."))
end

function get_mass_matrix(gp::GridapProblem, x, p, dΩ = gp.dΩ)
    return _get_mass_matrix(gp.mass_type, gp, x, p, dΩ)
end

_get_mass_matrix(::MassDefaut, gp::GridapProblem, x, p, dΩ) =
    _get_mass_matrix(MassDefaut(), gp, dΩ)

_get_mass_matrix(::ConstMass, gp::GridapProblem, x, p, dΩ) =
    _assemble_mass_matrix(gp, gp.mass)

function _get_mass_matrix(::SDMass, gp::GridapProblem, x, p, dΩ)
    uh = FEFunction(gp.U, x)
    return _assemble_mass_matrix(gp, (u, v) -> gp.mass(uh, p, u, v))
end
################################################################################
# structure to help casting the functional in a way that BifurcationKit can use
"""
    GridapBifProblem(res, u0, parms, V, U, dΩ, lens; kwargs...)

Construct a bifurcation problem which encodes a system of PDEs discretized with
[`Gridap`](https://github.com/gridap/Gridap.jl). The returned object is a subtype of
`BifurcationKit.AbstractDAEBifProblem`, so that a (possibly singular) mass matrix can be
used by the stability / eigenvalue solvers, *e.g.* for Hopf bifurcations of
differential-algebraic systems. It can be passed to `BifurcationKit.solve`,
`BifurcationKit.continuation`, `BifurcationKit.get_normal_form`, etc.

# Arguments
- `res(u, p, v)`: residual of the (semi-)discretized problem for the state `u`, the
  parameters `p` and the test function `v`. This is a Gridap weak form returning a
  `DomainContribution`.
- `u0`: initial guess, given either as a Gridap `FEFunction` (its free dof values are
  extracted) or directly as the vector of free dof values.
- `parms`: the set of parameters, typically a `NamedTuple` (*e.g.* `(λ = 1.0,)`) or a
  `Vector`.
- `V`: the `TestFESpace`.
- `U`: the `TrialFESpace` (it carries the Dirichlet data).
- `dΩ`: the `Measure` used to assemble the residual. It is stored and reused to assemble
  the mass matrix (see [`get_mass_matrix`](@ref)).
- `lens`: an `Accessors` optic selecting the continuation parameter inside `parms`, *e.g.*
  `(@optic _.λ)`, or `(@optic _[1])` for a vector of parameters.

# Keyword arguments
- `jac(u, p, du, v)`: analytical jacobian. If `nothing`, it is computed by finite differences.
- `d2res(u, p, du1, du2, v)`, `d3res(u, p, du1, du2, du3, v)`: second and third derivatives,
  required for automatic branch switching with a non-simple kernel.
- `mass`: the mass operator. It can be
  * `nothing`: the default L² mass `∫(u⋅v)*dΩ` (`MassDefaut`);
  * `mass(u, v)`: a constant, state-independent integrand (`ConstMass`), *e.g.*
    `(u,v) -> ∫(u⋅v)*dΩ`;
  * `mass(u, p, du, v)`: a state-dependent operator `M(x, p)` (`SDMass`), where
    `u` is the current state, `du` the trial function and `v` is the test function.
  See [`get_mass_matrix`](@ref).
- `record_from_solution`, `plot_solution`, `delta`: see the `BifurcationKit`
  documentation.
- `R01`, `R02`, `R11`: parameter-derivative operators. In addition to the
  `BifurcationKit` sentinels (`FiniteDifferences()`, `AutoDiff()`, `nothing`) and to
  an already assembled closure `(x, p) -> Vector`, they accept **analytic** Gridap
  weak forms, which are assembled automatically with the same convention as `res`
  and `jac`:
  * `R01(u, p, v)`, `R02(u, p, v)`: residual-like forms (return a
    `DomainContribution`), assembled into the free dof vector `∂_p F` (resp.
    `∂²_p F`);
  * `R11(u, p, du, v)`: jacobian-like form, assembled into the matrix `∂_p J` and
    applied to the direction `dx`.
  In these closures, `p` is the **full parameter `NamedTuple`** (`p.λ`, ...), *not*
  the scalar value of the continuation parameter.

# Extended methods
- [`get_mass_matrix`](@ref) assembles the mass matrix associated to `mass`.
- `BifurcationKit.is_mass_matrix_constant` and `BifurcationKit.getmassmatrix` are
  specialized so that the mass matrix is used by the DAE eigensolvers.

# See also
- [`GridapProblem`](@ref), [`get_mass_matrix`](@ref)
- `BifurcationKit.AbstractDAEBifProblem`, `BifurcationKit.solve`
"""
struct GridapBifProblem{Tfe, Tu, Tp, Tl, Tplot, Trec, Tδ, Tjet} <: BifurcationKit.AbstractDAEBifProblem
    "gridap problem"
    probFE::Tfe
    "Initial guess"
    u0::Tu
    "parameters"
    params::Tp
    "Typically a `Accessors.PropertyLens`. It specifies which parameter axis among `params` is used for continuation. For example, if `par = (α = 1.0, β = 1)`, we can perform continuation w.r.t. `α` by using `lens = (@optic _.α)`. If you have an array `par = [ 1.0, 2.0]` and want to perform continuation w.r.t. the first variable, you can use `lens = (@optic _[1])`. For more information, we refer to `Accessors.jl`."
    lens::Tl
    "user function to plot solutions during continuation. Signature: `plotSolution(x, p; kwargs...)`"
    plotSolution::Tplot
    "`record_from_solution = (x, p) -> norm(x)` function used record a few indicators about the solution. It could be `norm` or `(x, p) -> x[1]`. This is also useful when saving several huge vectors is not possible for memory reasons (for example on GPU...). This function can return pretty much everything but you should keep it small. For example, you can do `(x, p) -> (x1 = x[1], x2 = x[2], nrm = norm(x))` or simply `(x, p) -> (sum(x), 1)`. This will be stored in `contres.branch` (see below). Finally, the first component is used to plot in the continuation curve."
    recordFromSolution::Trec
    "used internally to compute derivatives (with finite differences) w.r.t the parameter `p`."
    δ::Tδ
    "Taylor jet w.r.t. parameters."
    jet::Tjet
end

# constructors (see docstring of the `GridapBifProblem` type above)
function GridapBifProblem(res, u0, parms, V, U, dΩ, lens;
                autodiff = false,
                jac = nothing,
                d2res = nothing,
                d3res = nothing,
                record_from_solution = BK.record_sol_default,
                plot_solution = BK.plot_default,
                R01 = BK.FiniteDifferences(),
                R02 = BK.FiniteDifferences(),
                R11 = BK.FiniteDifferences(),
                delta = BK._getprecision(Gridap.get_free_dof_values(u0)),
                mass = nothing,
                kwargs_jet...)
    jacFE =  autodiff ? nothing : jac
    probFE = GridapProblem(res, jacFE, d2res, d3res, V, U, nothing, dΩ, mass, _mass_type(mass))
    # analytic parameter-derivative forms (if provided) are assembled on the fly,
    # with the same convention as `res`/`jac`; BifurcationKit sentinels and
    # already assembled `(x, p) -> Vector` closures are forwarded untouched.
    R01jet = _is_residual_form(R01) ? ((x, p)     -> _residual_from_form(probFE, R01, x, p)) : R01
    R02jet = _is_residual_form(R02) ? ((x, p)     -> _residual_from_form(probFE, R02, x, p)) : R02
    R11jet = _is_jac_form(R11)      ? ((x, p, dx) -> _apply_from_biform(probFE, R11, x, p, dx)) : R11
    # type unstable but simplifies the types a lot
    jet = BK.Jet(; δ = delta, R01 = R01jet, R02 = R02jet, R11 = R11jet, kwargs_jet...)
    return GridapBifProblem(probFE, Gridap.get_free_dof_values(u0), parms, lens, plot_solution, record_from_solution, delta, jet)
end

BifurcationKit._getvectortype(gp::GridapProblem) = Gridap.FESpaces.get_vector_type(gp.U)
BifurcationKit._getvectortype(pb::GridapBifProblem) = BifurcationKit._getvectortype(pb.probFE)
BifurcationKit.isinplace(::GridapBifProblem) = false
BifurcationKit.residual(pb::GridapBifProblem, u, p) = residual(pb.probFE, u, p)
BifurcationKit.jacobian(pb::GridapBifProblem, u, p) = jacobian(pb.probFE, u, p)
BifurcationKit.dF(pb::GridapBifProblem, u, p, dx) = BifurcationKit.apply(BifurcationKit.jacobian(pb, u, p), dx)

BifurcationKit.d2F(pb::GridapBifProblem, u, p, dx1::AbstractArray{<:Real}, dx2::AbstractArray{<:Real}) = pb.probFE(u, p, dx1, dx2)
function BifurcationKit.d2F(pb::GridapBifProblem, x, p, dx1, dx2)
    probFE = pb.probFE
    dx1r = real.(dx1); dx2r = real.(dx2)
    dx1i = imag.(dx1); dx2i = imag.(dx2)
    return probFE(x, p, dx1r, dx2r) .- 
           probFE(x, p, dx1i, dx2i) .+ 
           im .* (probFE(x, p, dx1r, dx2i) .+ 
                  probFE(x, p, dx1i, dx2r))
end

BifurcationKit.d3F(pb::GridapBifProblem, u, p, dx1, dx2, dx3) = pb.probFE(u, p, dx1, dx2, dx3)
BifurcationKit.is_symmetric(::GridapBifProblem) = false
BifurcationKit.has_adjoint(::GridapBifProblem) = false
BifurcationKit.getdelta(pb::GridapBifProblem) = pb.δ
BifurcationKit.save_solution(::GridapBifProblem, x, p) = x
BifurcationKit.has_adjoint_MF(::GridapBifProblem) = false # TODO improve this using AD
BifurcationKit.update!(::GridapBifProblem, args...) = true
BifurcationKit.residual!(prob::GridapBifProblem, out, x, p) = out .= BK.residual(prob, x, p)

function Base.show(io::IO, prob::GridapBifProblem; prefix = "")
    gp = prob.probFE
    print(io, prefix, "┌─ Gridap Bifurcation Problem with uType ")
    printstyled(io, typeof(prob.u0), color = :cyan, bold = true)
    print(io, "\n", prefix, "├─ Inplace    : ")
    printstyled(io, false, color = :cyan, bold = true)
    print(io, "\n", prefix, "├─ Dimension  : ")
    printstyled(io, length(prob.u0), color = :cyan, bold = true)
    print(io, "\n", prefix, "├─ Parameter  : ")
    printstyled(io, BK.get_lens_symbol(BK.getlens(prob)), color = :cyan, bold = true)
    print(io, " = ", BK.getparam(prob), "\n")
    println(io, prefix, "├─ Mass       : ", nameof(typeof(gp.mass_type)))
    println(io, prefix, "├─ Jacobian   : ", gp.jac === nothing ? "finite differences" : "analytic")
    println(io, prefix, "└─ Spaces:    test  = ", nameof(typeof(gp.V)),
          ",\n              trial = ", nameof(typeof(gp.U)))
end

get_mass_matrix(prob::GridapBifProblem) = get_mass_matrix(prob.probFE)
get_mass_matrix(prob::GridapBifProblem, x, p) = get_mass_matrix(prob.probFE, x, p)
BK.has_hessian(::GridapBifProblem) = true

BK.R01(prob::GridapBifProblem, x, p) = BK.R01(BK.has_R01_trait(prob.jet), prob, x, p)
BK.R01(::BK.TraitUserPassed, prob::GridapBifProblem, x, p) = prob.jet.R01(x, p)
BK.R02(prob::GridapBifProblem, x, p) = BK.R02(BK.has_R02_trait(prob.jet), prob, x, p)
BK.R02(::BK.TraitUserPassed, prob::GridapBifProblem, x, p) = prob.jet.R02(x, p)

BK.R11(prob::GridapBifProblem, x, p, dx) = BK.R11(BK.has_R11_trait(prob.jet), prob, x, p, dx)
BK.R11(::BK.TraitUserPassed, prob::GridapBifProblem, x, p, dx) = prob.jet.R11(x, p, dx)

BK.is_mass_matrix_constant(prob::GridapBifProblem) = _is_constant_mass(prob.probFE.mass_type)
BK.getmassmatrix(prob::GridapBifProblem, x, p) = get_mass_matrix(prob.probFE, x, p)