import BifurcationKit: _getvectortype

# ─────────────────────────────────────────────────────────────────────────────
# Mass operator kinds, form detection and ddl helpers
# ─────────────────────────────────────────────────────────────────────────────

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

# A user form is recognised by its arguments: `methods` counts the callable itself,
# so a 3/4/5-argument form is reported as `nargs == 4/5/6`.
_has_arity(f, n) = f isa Function && any(m -> m.nargs == n, methods(f))

# Analytic parameter-derivative forms: `R01(u, p, v)` / `R02(u, p, v)` are
# residual-like, `R11(u, p, du, v)` is jacobian-like. Anything else (BK
# sentinels `FiniteDifferences()`/`AutoDiff()`, `nothing`, or an already
# assembled `(x, p) -> Vector` closure) is forwarded untouched.
_is_residual_form(f)  = _has_arity(f, 4)
_is_jac_form(f)       = _has_arity(f, 5)

# Analytic state-derivative of a state-dependent mass: a `dM(u, p, du1, du2, v)`
# trilinear form has arity 6.
_is_mass_hess_form(f) = _has_arity(f, 6)

# Build an FE function whose (constant) Dirichlet values share the dual eltype
# of the free values. Gridap otherwise rejects mixing dual free values with
# Float64 Dirichlet values.
function _dual_fe_function(f::MultiFieldFESpace, z)
    dir = [eltype(z).(Gridap.FESpaces.get_dirichlet_dof_values(sp)) for sp in f.spaces]
    return FEFunction(f, z, dir)
end
function _dual_fe_function(f::FESpace, z)
    return FEFunction(f, z, eltype(z).(Gridap.FESpaces.get_dirichlet_dof_values(f)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameter-agnostic FE problem
# ─────────────────────────────────────────────────────────────────────────────

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
  direction `du`. If `nothing`, the jacobian is built by `Gridap` from `res`
  (automatic differentiation) when the corresponding `FEOperator` is built.
- `d2res(u, p, du1, du2, v)`: second derivative (Hessian) of `res`. May also be
  a `BK` sentinel (`BK.FiniteDifferences()`, `BK.AutoDiff()`), in which case it
  is reconstructed accordingly.
- `d3res(u, p, du1, du2, du3, v)`: third derivative (Tressian) of `res`. May
  also be a `BK` sentinel (`BK.FiniteDifferences()`, `BK.AutoDiff()`).
- `V`: the `TestFESpace`.
- `U`: the `TrialFESpace` (it carries the Dirichlet data).
- `ls`: reserved linear solver, or `nothing` (currently unused).
- `dΩ`: the `Measure` used to assemble the residual; it is stored so that the
  mass matrix can be assembled later (see [`get_mass_matrix`](@ref)).
- `mass`: the mass operator, wrapped in a `BK.MassFunction` (or `nothing` for the
  default L² mass). The `MassFunction` carries the mass together with the
  derivatives required by the minimally augmented formulations:
  * `mass.M`: the constant integrand `mass(u, v)` or the state/parameter-dependent
    form `mass(u, p, du, v)` (or `nothing`);
  * `mass.Mᵗ`: an optional dedicated adjoint of the mass (unused here, `nothing`);
  * `mass.R01`: the operator `(x, p, v, w) -> ∂_p ⟨w, M(x,p) v⟩` (default
    `BK.AutoDiff()`), see `BK.R01_mass_matrix`;
  * `mass.∇x`: the operator `(x, p, v, w) -> ∇_x ⟨w, M(x,p) v⟩` (default
    `BK.AutoDiff()`), see `BK.∇_x_mass_matrix`.
  See [`get_mass_matrix`](@ref) and the *Mass matrix and derivatives* section
  below for the accepted signatures.
- `mass_type`: kind of mass operator, one of `MassDefaut()` (default L² mass),
  `ConstMass()` (state-independent bilinear integrand) or `SDMass()`
  (state-dependent `mass(u, p, du, v)`). It is inferred from `mass` by
  `_mass_type`.

# Mass matrix and derivatives
The mass matrix is assembled on the free dofs by
[`get_mass_matrix`](@ref) `(gp[, x, p])`; the state `x` is required for a
state-dependent mass (`SDMass`), while a constant mass (`MassDefaut`,
`ConstMass`) can be assembled without it. The kind is reported by `mass_type`
and, at the level of a [`GridapBifProblem`](@ref), by
`BK.is_mass_matrix_constant`.

When the mass depends on the state or on the parameters, the minimally augmented
Hopf formulation (`BK.newton_hopf`, `BK.continuation_hopf`) needs its
derivatives, exposed through
- `BK.R01_mass_matrix(gp, x, p, v, w)`: returns ``∂_p ⟨w, M(x,p) v⟩``,
  differentiating along the continuation lens;
- `BK.∇_x_mass_matrix(gp, x, p, v, w)`: returns ``∇_x ⟨w, M(x,p) v⟩``.

Both dispatch on the corresponding fields of `gp.mass`:
* `mass.R01` / `mass.∇x` given as a user closure `(x, p, v, w) -> scalar` /
  `(x, p, v, w) -> vector` are used as is;
* the `BK.AutoDiff()` and `BK.FiniteDifferences()` sentinels are evaluated by
  ForwardDiff (resp. central finite differences) on `get_mass_matrix(gp, x, p)`.

They are ignored for a constant mass, and are assembled with dual-valued free
dofs when needed (Dirichlet values are lifted to the dual type accordingly).

# Callable / extended methods
- `(gp::GridapProblem)(u, p, du1, du2)` returns the vector associated to
  `d2res(u, p, du1, du2, v)` (used by `BK.d2F`).
- `(gp::GridapProblem)(u, p, du1, du2, du3)` returns the vector associated to
  `d3res(u, p, du1, du2, du3, v)` (used by `BK.d3F`).
- `residual` and `jacobian` evaluate `res` and `jac` at a given state and
  parameter value.
- [`get_mass_matrix`](@ref) assembles the mass matrix, `BK.R01_mass_matrix` and
  `BK.∇_x_mass_matrix` its parameter and state derivatives (see above).

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

# short description of a derivative operator (jet entries `R01`/`R02`/`R11`,
# mass derivatives `R01M`/`∇xM`)
_deriv_op_name(::BK.AutoDiff) = "AutoDiff"
_deriv_op_name(::BK.FiniteDifferences) = "FiniteDifferences"
_deriv_op_name(::Nothing) = "missing"
_deriv_op_name(f) = "user"

# effective backend reported by `show`: `nothing` falls back to Gridap's
# automatic differentiation for `jac` and to finite differences for
# `d2res`/`d3res` (see `_d2F_fd`/`_d3F_fd`); user forms are reported as `user`.
_jac_name(jac) = jac === nothing ? "AutoDiff" : _deriv_op_name(jac)
_deriv_name(d) = d === nothing ? "FiniteDifferences" : _deriv_op_name(d)

# the mass derivatives only matter for a (state/parameter) dependent mass, in
# which case they are reported as `(R01=…, ∇x=…)`.
function _mass_deriv_str(gp::GridapProblem)
    _is_constant_mass(gp.mass_type) && return ""
    mf = gp.mass
    mf isa BK.MassFunction || return ""
    return string(" (R01=", _deriv_op_name(mf.R01), ", ∇x=", _deriv_op_name(mf.∇x), ")")
end

function Base.show(io::IO, gp::GridapProblem; prefix = "")
    print(io, prefix, "Gridap FE problem\n")
    print(io, prefix, "├─ Trial/Test : ", nameof(typeof(gp.U)), " / ", nameof(typeof(gp.V)), "\n")
    print(io, prefix, "├─ Free dofs  : ", Gridap.FESpaces.num_free_dofs(gp.U), "\n")
    print(io, prefix, "├─ Mass       : ", nameof(typeof(gp.mass_type)), _mass_deriv_str(gp), "\n")
    print(io, prefix, "└─ Derivatives: jac=", _jac_name(gp.jac),
          ", d2res=", _deriv_name(gp.d2res), ", d3res=", _deriv_name(gp.d3res))
end

# ─── operator rebuilt for a given parameter value
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

function jacobian!(A, gp::GridapProblem, u, p)
    op = op_from_param(gp, p)
    algop = Gridap.FESpaces.get_algebraic_operator(op)
    return Gridap.FESpaces.jacobian!(A, algop, u)
end

# ─── assembly of the analytic parameter-derivative forms
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

# ─── derivative flavours: analytic form, finite differences, ForwardDiff
# FE function associated to a *direction* `d` (a free dof vector). The Hessian
# and Tressian are derivatives w.r.t. the free dofs, so the direction must carry
# *homogeneous* Dirichlet data: `FEFunction(U, d)` would instead inject the
# (generally non-zero) Dirichlet data of `U`.
_homogeneous_fe(gp::GridapProblem, d) = FEFunction(gp.U, d, Gridap.FESpaces.zero_dirichlet_values(gp.U))

# assemble the gradient `∇_x ⟨w, M(x,p) v⟩` from the analytic trilinear mass
# derivative form `dM(u, p, du1, du2, v)` (used for `∇xM`). The direction `du2`
# is fixed to `v` (the mass direction), `du1` is the state-derivative trial and
# the test index is contracted with `w`. This mirrors
# `_∇_x_mass_matrix(::BK.AutoDiff, …)` with its real / imaginary split.
function _mass_grad_from_form(U, V, form, x, p, v, w)
    uh = _dual_fe_function(U, x)
    function _grad(a, b)
        ah = FEFunction(U, a, Gridap.FESpaces.zero_dirichlet_values(U))
        A  = Gridap.FESpaces.assemble_matrix((du1, test) -> form(uh, p, du1, ah, test), U, V)
        return transpose(A) * b
    end
    vr = real(v); vi = imag(v)
    wr = real(w); wi = imag(w)
    return complex.(_grad(vr, wr) .+ _grad(vi, wi),
                    _grad(vi, wr) .- _grad(vr, wi))
end

# second derivative from the analytic form `d2res`
function (gp::GridapProblem)(u, p, du1, du2)
    du1h = _homogeneous_fe(gp, du1)
    du2h = _homogeneous_fe(gp, du2)
    a(u, v) = gp.d2res(u, p, du1h, du2h, v)
    feop = FEOperator(a, gp.U, gp.V)
    alop = Gridap.FESpaces.get_algebraic_operator(feop)
    Gridap.FESpaces.residual(alop, u)
end

# third derivative from the analytic form `d3res`
function (gp::GridapProblem)(u, p, du1, du2, du3)
    du1h = _homogeneous_fe(gp, du1)
    du2h = _homogeneous_fe(gp, du2)
    du3h = _homogeneous_fe(gp, du3)
    a(u, v) = gp.d3res(u, p, du1h, du2h, du3h, v)
    feop = FEOperator(a, gp.U, gp.V)
    alop = Gridap.FESpaces.get_algebraic_operator(feop)
    Gridap.FESpaces.residual(alop, u)
end

# second/third derivative by finite differences (`d2res`/`d3res` absent)
_d2F_fd(gp::GridapProblem, u, p, du1, du2) = jvp(central_fdm(3, 1), z -> jacobian(gp, z, p) * du1, (u, du2))
_d3F_fd(gp::GridapProblem, u, p, du1, du2, du3) = jvp(central_fdm(3, 1), z -> gp(z, p, du1, du2), (u, du3))

# Assemble the residual at dual-typed dof values `z` into a free-dof vector
# allocated with the dual eltype (the standard assembler would allocate a
# Float vector).
function _residual_dual(gp::GridapProblem, p, z)
    uh = _dual_fe_function(gp.U, z)
    v  = get_fe_basis(gp.V)
    data = Gridap.FESpaces.collect_cell_vector(gp.V, gp.res(uh, p, v))
    b = similar(Gridap.get_free_dof_values(zero(gp.V)), eltype(z))
    fill!(b, zero(eltype(z)))
    Gridap.FESpaces.assemble_vector!(b, SparseMatrixAssembler(gp.U, gp.V), data)
    return b
end

# second/third derivative by ForwardDiff on the assembled residual (never on `jacobian`)
function _d2F_ad(gp::GridapProblem, u::AbstractArray{𝒯1}, p, du1::AbstractArray{𝒯2}, du2) where {𝒯1, 𝒯2}
    𝒯 = promote_type(𝒯1, 𝒯2)
    return ForwardDiff.derivative(
        ε2 -> ForwardDiff.derivative(
            ε1 -> _residual_dual(gp, p, u .+ ε1 .* du1 .+ ε2 .* du2), zero(𝒯)), zero(𝒯))
end
function _d3F_ad(gp::GridapProblem, u::AbstractArray{𝒯1}, p, du1::AbstractArray{𝒯2}, du2, du3) where {𝒯1, 𝒯2}
    𝒯 = promote_type(𝒯1, 𝒯2)
    return ForwardDiff.derivative(
        ε3 -> ForwardDiff.derivative(
            ε2 -> ForwardDiff.derivative(
                ε1 -> _residual_dual(gp, p, u .+ ε1 .* du1 .+ ε2 .* du2 .+ ε3 .* du3), zero(𝒯)), zero(𝒯)), zero(𝒯))
end

(gp::GridapProblem{Tres, Tjac, Nothing})(u, p, du1, du2) where {Tres, Tjac} = _d2F_fd(gp, u, p, du1, du2)
(gp::GridapProblem{Tres, Tjac, BK.FiniteDifferences})(u, p, du1, du2) where {Tres, Tjac} = _d2F_fd(gp, u, p, du1, du2)
(gp::GridapProblem{Tres, Tjac, BK.AutoDiff})(u, p, du1, du2) where {Tres, Tjac} = _d2F_ad(gp, u, p, du1, du2)

(gp::GridapProblem{Tres, Tjac, Td2res, Nothing})(u, p, du1, du2, du3) where {Tres, Tjac, Td2res} = _d3F_fd(gp, u, p, du1, du2, du3)
(gp::GridapProblem{Tres, Tjac, Td2res, BK.FiniteDifferences})(u, p, du1, du2, du3) where {Tres, Tjac, Td2res} = _d3F_fd(gp, u, p, du1, du2, du3)
(gp::GridapProblem{Tres, Tjac, Td2res, BK.AutoDiff})(u, p, du1, du2, du3) where {Tres, Tjac, Td2res} = _d3F_ad(gp, u, p, du1, du2, du3)

# ─────────────────────────────────────────────────────────────────────────────
# Mass matrix assembly
# ─────────────────────────────────────────────────────────────────────────────

mass_default(u,v,dΩ) = ∫(u⋅v) * dΩ

# scalar entry type of the (possibly block) cell matrices, e.g.
# `Matrix{Float64} -> Float64` and `ArrayBlock{Matrix{Dual},2} -> Dual`.
_scalar_eltype(::Type{T}) where {T} = (U = eltype(T); U === T ? T : _scalar_eltype(U))

# shared assembly of a (bi)linear mass integrand over the free dofs
function _assemble_mass_matrix(gp::GridapProblem, integrand)
    (;U, V) = gp
    v = get_fe_basis(V)
    u = get_trial_fe_basis(U)
    assemblytuple = Gridap.FESpaces.collect_cell_matrix(U,V,integrand(u,v))
    cell_matrix_MM   = collect(assemblytuple[1][1]) # This result is no longer a LazyArray
    newassemblytuple = ([cell_matrix_MM], assemblytuple[2], assemblytuple[3])
    # the scalar type is read off the cell matrices so that dual-valued states
    # (ForwardDiff derivatives of a state-dependent mass) are assembled in the
    # dual type instead of the default `Float64`.
    T = _scalar_eltype(eltype(cell_matrix_MM))
    a = SparseMatrixAssembler(Gridap.Algebra.SparseMatrixCSC{T,Int}, Vector{T}, U, V)
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
get_mass_matrix(gp::GridapProblem, dΩ = gp.dΩ) = _get_mass_matrix(gp.mass_type, gp, dΩ)
_get_mass_matrix(::MassDefaut, gp::GridapProblem, dΩ) = _assemble_mass_matrix(gp, (u, v) -> mass_default(u, v, dΩ))
_get_mass_matrix(::ConstMass, gp::GridapProblem, dΩ) = _assemble_mass_matrix(gp, gp.mass.M)

get_mass_matrix(gp::GridapProblem, x, p, dΩ = gp.dΩ) = _get_mass_matrix(gp.mass_type, gp, x, p, dΩ)
_get_mass_matrix(::MassDefaut, gp::GridapProblem, x, p, dΩ) = _get_mass_matrix(MassDefaut(), gp, dΩ)
_get_mass_matrix(::ConstMass, gp::GridapProblem, x, p, dΩ) = _assemble_mass_matrix(gp, gp.mass.M)

function _get_mass_matrix(::SDMass, ::GridapProblem, dΩ)
    throw(ArgumentError("The mass operator passed to `GridapBifProblem` is state-dependent; call `get_mass_matrix(gp, x, p)` with the state `x` and parameters `p`."))
end

function _get_mass_matrix(::SDMass, gp::GridapProblem, x, p, dΩ)
    # `_dual_fe_function` keeps the dual eltype of `x` on the Dirichlet values,
    # which Gridap otherwise rejects (needed by the AutoDiff mass derivatives).
    uh = _dual_fe_function(gp.U, x)
    return _assemble_mass_matrix(gp, (u, v) -> gp.mass.M(uh, p, u, v))
end

# ─────────────────────────────────────────────────────────────────────────────
# Mass derivatives, expected by the minimally augmented Hopf formulation of BK.
# BK only defines `R01_mass_matrix` / `∇_x_mass_matrix` for `DAEMassBifProblem`.
# A `GridapBifProblem` stores its mass in `probFE.mass` (a `BK.MassFunction`);
# we reproduce here the BK sentinel dispatch (user form / finite differences /
# ForwardDiff) on top of the Gridap-aware `getmassmatrix`.
# `R01_mass_matrix(prob, x, p, v, w)` returns `∂_p ⟨w, M(x, p) v⟩` and
# `∇_x_mass_matrix(prob, x, p, v, w)` returns `∇_x ⟨w, M(x, p) v⟩`.
# The `GridapBifProblem` entry points are defined with the BK interface below.
# ─────────────────────────────────────────────────────────────────────────────

# user provided form `(x, p, v, w) -> scalar`
_R01_mass_matrix(R01, prob, x, p, v, w) = R01(x, p, v, w)

function _R01_mass_matrix(::BK.FiniteDifferences, prob, x, p, v, w)
    lens = BK.getlens(prob)
    p0 = BK._get(p, lens)
    ϵ = BK.getdelta(prob)
    M₊ = BK.getmassmatrix(prob, x, BK.set(p, lens, p0 + ϵ))
    M₋ = BK.getmassmatrix(prob, x, BK.set(p, lens, p0 - ϵ))
    return (BK.dot_with_mass(w, M₊, v) - BK.dot_with_mass(w, M₋, v)) / (2ϵ)
end

function _R01_mass_matrix(::BK.AutoDiff, prob, x, p, v, w)
    lens = BK.getlens(prob)
    p0 = BK._get(p, lens)
    return ForwardDiff.derivative(z -> BK.dot_with_mass(w, BK.getmassmatrix(prob, x, BK.set(p, lens, z)), v), p0)
end

# user provided form `(x, p, v, w) -> vector`
_∇_x_mass_matrix(∇x, prob, x, p, v, w) = ∇x(x, p, v, w)

function _∇_x_mass_matrix(::BK.AutoDiff, prob, x, p, v, w)
    # real / imaginary split to stay within ForwardDiff's real arithmetic
    vr = real(v); vi = imag(v)
    wr = real(w); wi = imag(w)
    quad(a, b, z) = BK.dot_with_mass(b, BK.getmassmatrix(prob, z, p), a)
    gre = ForwardDiff.gradient(z -> quad(vr, wr, z) + quad(vi, wi, z), x)
    gim = ForwardDiff.gradient(z -> quad(vi, wr, z) - quad(vr, wi, z), x)
    return complex.(gre, gim)
end

# ─────────────────────────────────────────────────────────────────────────────
# BK-facing bifurcation problem
# ─────────────────────────────────────────────────────────────────────────────

"""
    GridapBifProblem(res, u0, parms, V, U, dΩ, lens; kwargs...)

Encode a system of PDEs discretized with [`Gridap`](https://github.com/gridap/Gridap.jl)
as a bifurcation problem of the (semi-discretized) differential-algebraic form

```math
M(x,p)\\, \\dot x = F(x,p),
```

where the vector field ``F`` is given by the weak form `res` and the (possibly
singular) mass matrix ``M`` is assembled from the `mass` keyword (the L² mass by
default, see [`get_mass_matrix`](@ref)). The returned object is a subtype of
`BK.AbstractDAEBifProblem` and can be passed to `BK.solve`, `BK.continuation`,
`BK.newton_hopf`, `BK.get_normal_form`, etc. The mass matrix is used by the DAE
eigensolvers to assess stability and to detect Hopf points through the generalized
eigenproblem ``dF(x,p)\\,\\phi = \\lambda\\, M(x,p)\\,\\phi`` (see the
[Mass matrix](@ref mass-matrix) section).

# Arguments
- `res(u, p, v)`: residual of the (semi-)discretized problem, a weak form
  returning a `DomainContribution`. `u` is a `FEFunction` on `U` built from the
  current state, `p` the parameters and `v` a test function.
- `u0`: initial guess, either a `FEFunction` on `U` (its free dof values are
  extracted) or directly the vector of free dof values.
- `parms`: the parameters, typically a `NamedTuple` (*e.g.* `(λ = 1.0,)`) or a
  `Vector`.
- `V`: the `TestFESpace`.
- `U`: the `TrialFESpace` (it carries the Dirichlet data).
- `dΩ`: the `Measure` used to assemble the forms; stored and reused to assemble
  the mass matrix.
- `lens`: an `Accessors` optic selecting the continuation parameter in `parms`,
  *e.g.* `(@optic _.λ)`, or `(@optic _[1])` for a `Vector` of parameters.

# Keyword arguments
- `jacobian_type = BK.FullSparse()`: selects how `BK.jacobian(prob, x, p)` is
  evaluated. Accepted values:
  * `BK.FullSparse()` (default): a fresh sparse matrix is assembled at each call;
  * `BK.FullSparseInplace()`: a sparse matrix is preallocated when the problem is
    built and updated in place at each call (faster, but the sparsity pattern of
    the vector field must be constant);
  * `BK.MatrixFree()`: `BK.jacobian` returns a closure `dx -> J(x,p)*dx` computed
    by ForwardDiff on the residual, without assembling the matrix. It must be
    used together with a matrix-free linear solver (*e.g.* `GMRESIterativeSolvers`).
- `jac(u, p, du, v)`: analytic jacobian. If `nothing` (default), it is built by
  Gridap from `res` (automatic differentiation).
- `autodiff = false`: if `true`, `jac` is ignored and the jacobian is built by
  Gridap from `res`.
- `d2res`, `d3res`: second and third derivatives of the residual, used by
  `d2F`/`d3F` (*e.g.* for automatic branch switching with a non-simple kernel).
  Accepted values:
  * an analytic weak form `d2res(u, p, du1, du2, v)` /
    `d3res(u, p, du1, du2, du3, v)`, assembled like `res` (the directions `dui`
    are homogeneous, *i.e.* with zero Dirichlet data);
  * `BK.FiniteDifferences()` or `nothing` (default): finite differences;
  * `BK.AutoDiff()`: ForwardDiff applied to the assembled residual.
- `mass`: the mass operator. Accepted values:
  * `nothing` (default): the L² mass `∫(u⋅v)*dΩ`;
  * `mass(u, v)`: a constant, state-independent integrand;
  * `mass(u, p, du, v)`: a state- and/or parameter-dependent integrand
    `M(x, p)`, where `u` is the current state, `du` the trial function and `v`
    the test function.
  See [`get_mass_matrix`](@ref).
- `R01M`, `∇xM`: derivatives of the mass, used by the minimally augmented Hopf
  formulation (`BK.newton_hopf`, `BK.continuation_hopf`) when `mass` depends on
  the parameters or on the state:
  * `R01M(x, p, v, w) -> scalar` computes ``∂_p\\langle w, M(x,p) v\\rangle``,
    differentiating along `lens` (or the sentinels `BK.AutoDiff()`,
    `BK.FiniteDifferences()`, `nothing`);
  * `∇xM` computes ``\\nabla_x\\langle w, M(x,p) v\\rangle``. It accepts
    `BK.AutoDiff()` (default), `BK.FiniteDifferences()`, `nothing`, an already
    assembled closure `(x, p, v, w) -> Vector`, or an **analytic** Gridap weak
    form `dM(u, p, du1, du2, v)`: the state-derivative direction of the mass
    bilinear form, with `du2` the mass direction and `v` the test function
    (assembled with homogeneous Dirichlet directions, like `d2res`).
  They are ignored for a constant mass.
- `R01`, `R02`, `R11`: parameter-derivative operators used for stability and
  normal-form computations. In addition to the `BK` sentinels
  (`BK.FiniteDifferences()`, `BK.AutoDiff()`, `nothing`) and to an already
  assembled closure `(x, p) -> Vector`, they accept **analytic** weak forms,
  assembled with the same convention as `res`/`jac`:
  * `R01(u, p, v)`, `R02(u, p, v)`: residual-like forms, assembled into the free
    dof vector ``∂_p F`` (resp. ``∂^2_p F``);
  * `R11(u, p, du, v)`: jacobian-like form, assembled into the matrix ``∂_p J``
    and applied to the direction `du`.
  In these closures `p` is the **full parameter set** (`p.λ`, ...), *not* the
  scalar value of the continuation parameter.
- `record_from_solution = BK.record_sol_default`, `plot_solution = BK.plot_default`:
  callbacks used during continuation, see the `BK` documentation.
  `record_from_solution(x, p)` returns a small `NamedTuple` of indicators stored
  along the branch; `plot_solution(x, p; kwargs...)` displays the solution.
- `delta`: finite-difference step used by `BK` for the parameter derivatives
  (default: `sqrt(eps)` of the state's scalar type).
- `kwargs_jet...`: further keyword arguments forwarded to `BK.Jet`.

# Fields
- `probFE::GridapProblem`: the parameter-agnostic FE problem, see
  [`GridapProblem`](@ref).
- `jacobianType`: the selected jacobian kind, one of `BK.FullSparse()`,
  `BK.FullSparseInplace()` or `BK.MatrixFree()`, see `jacobian_type`.
- `jacobian`: the preallocated sparse matrix used by `BK.FullSparseInplace()`
  (updated in place), or `nothing` for the other kinds.
- `u0`: the initial state (free dof vector).
- `params`, `lens`: the parameters and the continuation optic.
- `plotSolution`, `recordFromSolution`: the callbacks.
- `δ`: the finite-difference step.
- `jet::BK.Jet`: the parameter-derivative operators.

# Extended methods
- [`get_mass_matrix`](@ref) and `BK.getmassmatrix(prob, x, p)` assemble the mass
  matrix; `BK.is_mass_matrix_constant(prob)` tells whether it is constant.
- `residual(prob, x, p)`, `jacobian(prob, x, p)` and `BK.dF`, `BK.d2F`, `BK.d3F`
  evaluate the vector field and its derivatives. `BK.jacobian` dispatches on
  `jacobian_type`: it returns a sparse matrix for `BK.FullSparse()` /
  `BK.FullSparseInplace()` and a matrix-free closure `dx -> J(x,p)*dx` for
  `BK.MatrixFree()`.
- `BK.R01`, `BK.R02`, `BK.R11` evaluate the parameter derivatives.
- `BK.R01_mass_matrix` and `BK.∇_x_mass_matrix` provide the mass derivatives
  required by the minimally augmented Hopf formulation.
- The problem is not inplace, not symmetric and has no adjoint
  (`BK.isinplace`, `BK.is_symmetric`, `BK.has_adjoint`).

# See also
- [`GridapProblem`](@ref), [`get_mass_matrix`](@ref)
- `BK.AbstractDAEBifProblem`, `BK.solve`, `BK.continuation`, `BK.newton_hopf`
- [Mass matrix](@ref mass-matrix), [Simple Hopf point](@ref simple-hopf)

# Example
```julia
# -u'' + u + u³ = λ on (0,1), with a state-dependent mass
res(u, p, v)     = ∫(∇(u) ⋅ ∇(v) + u * v + u^3 * v - p.λ * v) * dΩ
jac(u, p, du, v) = ∫(∇(du) ⋅ ∇(v) + du * v + 3 * u^2 * du * v) * dΩ
mass(u, p, du, v) = ∫((1 + u^2) * du * v) * dΩ

prob = GridapBifProblem(res, u0, (λ = 1.0,), V, U, dΩ, (@optic _.λ);
                        jac = jac, mass = mass)

sol = BK.solve(prob, BK.Newton(), optn)
br  = BK.continuation(prob, BK.Natural(), opts)
```
"""
struct GridapBifProblem{Tfe, Tjac, Tjc, Tu, Tp, Tl, Tplot, Trec, Tδ, Tjet} <: BK.AbstractDAEBifProblem
    "gridap problem"
    probFE::Tfe
    "selected jacobian kind: `BK.FullSparse()`, `BK.FullSparseInplace()` or `BK.MatrixFree()`"
    jacobianType::Tjac
    "preallocated jacobian matrix for `BK.FullSparseInplace()`, `nothing` otherwise"
    jacobian::Tjc
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
                jacobian_type = BK.FullSparse(),
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
                R01M = BK.AutoDiff(),
                ∇xM  = BK.AutoDiff(),

                kwargs_jet...)
    jacFE = autodiff ? nothing : jac
    # an analytic mass state-derivative form `dM(u, p, du1, du2, v)` (passed as
    # `∇xM`) is assembled on the fly; sentinels (`BK.AutoDiff()`, …) and already
    # assembled `(x, p, v, w) -> Vector` closures are forwarded untouched.
    ∇xMjet = _is_mass_hess_form(∇xM) ?
        ((x, p, v, w) -> _mass_grad_from_form(U, V, ∇xM, x, p, v, w)) : ∇xM
    massfun = BK.MassFunction(mass, nothing, R01M, ∇xMjet)
    probFE = GridapProblem(res, jacFE, d2res, d3res, V, U, nothing, dΩ, massfun, _mass_type(mass))
    # analytic parameter-derivative forms (if provided) are assembled on the fly,
    # with the same convention as `res`/`jac`; BK sentinels and
    # already assembled `(x, p) -> Vector` closures are forwarded untouched.
    R01jet = _is_residual_form(R01) ? ((x, p)     -> _residual_from_form(probFE, R01, x, p)) : R01
    R02jet = _is_residual_form(R02) ? ((x, p)     -> _residual_from_form(probFE, R02, x, p)) : R02
    R11jet = _is_jac_form(R11)      ? ((x, p, dx) -> _apply_from_biform(probFE, R11, x, p, dx)) : R11
    # type unstable but simplifies the types a lot
    jet = BK.Jet(; δ = delta, R01 = R01jet, R02 = R02jet, R11 = R11jet, kwargs_jet...)
    x0 = Gridap.get_free_dof_values(u0)
    J = _init_jacobian(probFE, jacobian_type, x0, parms)
    return GridapBifProblem(probFE, jacobian_type, J, x0, parms, lens, plot_solution, record_from_solution, delta, jet)
end

# preallocate the jacobian only for the in-place flavour (a fresh matrix is
# assembled for `FullSparse`, and none for `MatrixFree`).
_init_jacobian(probFE::GridapProblem, ::BK.FullSparseInplace, x0, p) = jacobian(probFE, x0, p)
_init_jacobian(::GridapProblem, ::BK.AbstractJacobianType, x0, p) = nothing

# ─────────────────────────────────────────────────────────────────────────────
# BifurcationKit interface
# ─────────────────────────────────────────────────────────────────────────────

# ─── traits and spaces
BK._getvectortype(gp::GridapProblem) = Gridap.FESpaces.get_vector_type(gp.U)
BK._getvectortype(pb::GridapBifProblem) = BK._getvectortype(pb.probFE)
BK.isinplace(::GridapBifProblem) = false
BK.is_symmetric(::GridapBifProblem) = false
BK.has_adjoint(::GridapBifProblem) = false
BK.has_adjoint_MF(::GridapBifProblem) = false # TODO improve this using AD
BK.has_hessian(::GridapBifProblem) = true
BK.update!(::GridapBifProblem, args...) = true
BK.getdelta(pb::GridapBifProblem) = pb.δ

# ─── residual, jacobian and derivatives
BK.residual(pb::GridapBifProblem, u, p) = residual(pb.probFE, u, p)

# select the jacobian flavour chosen at construction (`jacobian_type`)
BK.jacobian(pb::GridapBifProblem, u, p) = _jacobian(pb, pb.jacobianType, u, p)

# `FullSparse`: a new sparse matrix is assembled at each call
_jacobian(pb::GridapBifProblem, ::BK.FullSparse, u, p) = jacobian(pb.probFE, u, p)

# `FullSparseInplace`: update the matrix preallocated at construction (the
# sparsity pattern is assumed constant) and return it
function _jacobian(pb::GridapBifProblem, ::BK.FullSparseInplace, u, p)
    jacobian!(pb.jacobian, pb.probFE, u, p)
    return pb.jacobian
end

# `MatrixFree`: return the jacobian-vector product as a closure, without ever
# assembling the matrix. It must be used with a matrix-free linear solver.
_jacobian(pb::GridapBifProblem, ::BK.MatrixFree, u, p) = dx -> _jvp_ad(pb, u, p, dx)

function _jacobian(pb::GridapBifProblem, jt, u, p)
    throw(ArgumentError("jacobian_type = $jt is not supported; use BK.FullSparse(), BK.FullSparseInplace() or BK.MatrixFree()."))
end

# matrix-free jacobian-vector product by ForwardDiff on the assembled residual
function _jvp_ad(pb::GridapBifProblem, u, p, dx)
    𝒯 = promote_type(eltype(u), eltype(dx))
    return ForwardDiff.derivative(ε -> _residual_dual(pb.probFE, p, u .+ ε .* dx), zero(𝒯))
end

BK.jacobian!(pb::GridapBifProblem, J, u, p) = jacobian!(J, pb.probFE, u, p)
BK.dF(pb::GridapBifProblem, u, p, dx) = BK.apply(BK.jacobian(pb, u, p), dx)

BK.d2F(pb::GridapBifProblem, u, p, dx1::AbstractArray{<:Real}, dx2::AbstractArray{<:Real}) = pb.probFE(u, p, dx1, dx2)
function BK.d2F(pb::GridapBifProblem, x, p, dx1, dx2)
    @error "******* d2F"
    probFE = pb.probFE
    dx1r = real.(dx1); dx2r = real.(dx2)
    dx1i = imag.(dx1); dx2i = imag.(dx2)
    return probFE(x, p, dx1r, dx2r) .- 
           probFE(x, p, dx1i, dx2i) .+ 
           im .* (probFE(x, p, dx1r, dx2i) .+ 
                  probFE(x, p, dx1i, dx2r))
end

BK.d3F(pb::GridapBifProblem, u, p, dx1, dx2, dx3) = pb.probFE(u, p, dx1, dx2, dx3)
BK.residual!(prob::GridapBifProblem, out, x, p) = out .= BK.residual(prob, x, p)
BK.save_solution(::GridapBifProblem, x, p) = x

# ─── parameter derivatives
BK.R01(prob::GridapBifProblem, x, p) = BK.R01(BK.has_R01_trait(prob.jet), prob, x, p)
BK.R01(::BK.TraitUserPassed, prob::GridapBifProblem, x, p) = prob.jet.R01(x, p)
BK.R02(prob::GridapBifProblem, x, p) = BK.R02(BK.has_R02_trait(prob.jet), prob, x, p)
BK.R02(::BK.TraitUserPassed, prob::GridapBifProblem, x, p) = prob.jet.R02(x, p)

BK.R11(prob::GridapBifProblem, x, p, dx) = BK.R11(BK.has_R11_trait(prob.jet), prob, x, p, dx)
BK.R11(::BK.TraitUserPassed, prob::GridapBifProblem, x, p, dx) = prob.jet.R11(x, p, dx)

# ─── mass matrix and derivatives
get_mass_matrix(prob::GridapBifProblem) = get_mass_matrix(prob.probFE)
get_mass_matrix(prob::GridapBifProblem, x, p) = get_mass_matrix(prob.probFE, x, p)
BK.is_mass_matrix_constant(prob::GridapBifProblem) = _is_constant_mass(prob.probFE.mass_type)
BK.getmassmatrix(prob::GridapBifProblem, x, p) = get_mass_matrix(prob.probFE, x, p)

# mass derivatives expected by the minimally augmented Hopf formulation of BK;
# the implementations live next to `get_mass_matrix` (see above).
BK.R01_mass_matrix(prob::GridapBifProblem, x, p, v, w) = _R01_mass_matrix(prob.probFE.mass.R01, prob, x, p, v, w)
BK.∇_x_mass_matrix(prob::GridapBifProblem, x, p, v, w) = _∇_x_mass_matrix(prob.probFE.mass.∇x, prob, x, p, v, w)

# ─── display
function Base.show(io::IO, prob::GridapBifProblem; prefix = "")
    gp = prob.probFE
    print(io, prefix, "┌─ Gridap Bifurcation Problem with uType ")
    printstyled(io, typeof(prob.u0), color = :cyan, bold = true)
    print(io, "\n", prefix, "├─ Inplace       : ")
    printstyled(io, false, color = :cyan, bold = true)
    print(io, "\n", prefix, "├─ Dimension     : ")
    printstyled(io, length(prob.u0), color = :cyan, bold = true)
    print(io, "\n", prefix, "├─ Parameter     : ")
    printstyled(io, BK.get_lens_symbol(BK.getlens(prob)), color = :cyan, bold = true)
    print(io, " = ", BK.getparam(prob), "\n")
    println(io, prefix, "├─ Mass          : ", nameof(typeof(gp.mass_type)), _mass_deriv_str(gp))
    println(io, prefix, "├─ Jacobian      : ", nameof(typeof(prob.jacobianType)), " (", _jac_name(gp.jac), ")")
    println(io, prefix, "├─ Derivatives   : jac=", _jac_name(gp.jac), ", d2res=", _deriv_name(gp.d2res), ", d3res=", _deriv_name(gp.d3res))
    println(io, prefix, "├─ R01, R02, R11 : ", (_deriv_op_name(prob.jet.R01), _deriv_op_name(prob.jet.R02), _deriv_op_name(prob.jet.R11)))
    println(io, prefix, "└─ Spaces:    test  = ", nameof(typeof(gp.V)),
          ",\n              trial = ", nameof(typeof(gp.U)))
end
