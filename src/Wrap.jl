
struct GridapProblem{Tres, Tjac, Td2res, Td3res, TV, TU, Tls, TOm, Tm}
    res::Tres        # res(u, p, v),                 residual
    jac::Tjac        # jac(u, p, du, v),             jacobian
    d2res::Td2res    # d2res(u, du1, du2, v)
    d3res::Td3res    # d3res(u, du1, du2, du3, v)
    V::TV
    U::TU
    ls::Tls
    dΩ::TOm
    mass::Tm
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

# # user_uh_to_cell_residual
# function (gp::GridapProblem)(::Val{:CellRes}, u, p)
#     op = op_from_param(gp, p)
#     algop = Gridap.FESpaces.get_algebraic_operator(op)
#     return Gridap.FESpaces.residual(algop, u)
# end

mass_default(u,v,dΩ) = ∫(u⋅v) * dΩ

"""
    get_mass_matrix(prob::GridapBifProblem)
    get_mass_matrix(prob::GridapProblem, dΩ = prob.dΩ)

Assemble the (sparse) mass matrix associated to the problem, on the free dofs of the
trial/test spaces (Dirichlet dofs are eliminated, consistently with the jacobian).

The integrand is the one passed to `GridapBifProblem(...; mass = ...)`, defaulting to the
L² mass `∫(u⋅v)*dΩ`. For an incompressible flow, a typical choice is the *velocity only*
mass `(u,p),(v,q) -> ∫(v⊙u)*dΩ`, which yields a singular mass matrix with a zero pressure
block, as required for the stability of the differential-algebraic system
``M\\dot z = F(z, p)``.
"""
function get_mass_matrix(prob::GridapProblem, dΩ = prob.dΩ)
    (;U, V) = prob
    mm = isnothing(prob.mass) ? (u,v) -> mass_default(u,v,dΩ) : prob.mass
    v = get_fe_basis(V)
    u = get_trial_fe_basis(U)
    assemblytuple = Gridap.FESpaces.collect_cell_matrix(U,V,mm(u,v))
    cell_matrix_MM   = collect(assemblytuple[1][1]) # This result is no longer a LazyArray
    newassemblytuple = ([cell_matrix_MM], assemblytuple[2], assemblytuple[3])
    a = SparseMatrixAssembler(U, V)
    L2MM = assemble_matrix(a, newassemblytuple)
    return L2MM
end
################################################################################
# structure to help casting the functional in a way that BifurcationKit can use
"""
    GridapBifProblem(res, u0, parms, V, U, dΩ, lens; jac = nothing, mass = nothing, kwargs...)

Construct a bifurcation problem which encodes a system of PDEs discretized with
[`Gridap`](https://github.com/gridap/Gridap.jl). It is a subtype of
`BifurcationKit.AbstractDAEBifProblem` so that a (possibly singular) mass matrix can be
used for the stability analysis, *e.g.* for Hopf bifurcations.

# Arguments
- `res(u, p, v)`: residual of the (semi-)discretized problem, where `p` are the parameters.
- `u0`: initial guess (a `Gridap` `FEFunction` or its free dof values).
- `parms`: the set of parameters.
- `V`: `TestFESpace`.
- `U`: `TrialFESpace`.
- `dΩ`: the `Measure` used to assemble the residual (stored for the mass matrix).
- `lens`: an `Accessors` lens selecting the continuation parameter in `parms`, *e.g.*
  `(@optic _.λ)`.

# Keyword arguments
- `jac(u, p, du, v)`: analytical jacobian. If `nothing`, it is computed by finite differences.
- `d2res(u, p, du1, du2, v)`, `d3res(u, p, du1, du2, du3, v)`: second and third derivatives,
  required for automatic branch switching with a non-simple kernel.
- `mass(u, v)`: integrand of the mass matrix, *e.g.* `(u,v) -> ∫(u⋅v)*dΩ`. If `nothing`,
  the default L² mass `∫(u⋅v)*dΩ` is used (see [`get_mass_matrix`](@ref)).
- `record_from_solution`, `plot_solution`, `R01`, `R02`, `R11`, `delta`: see the
  `BifurcationKit` documentation.

# Extended methods
- [`get_mass_matrix`](@ref) assembles the mass matrix associated to `mass`.
- `BifurcationKit.is_mass_matrix_constant` and `BifurcationKit.getmassmatrix` are
  specialized so that the mass matrix is used by the DAE eigensolvers.
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

import BifurcationKit: _getvectortype

BifurcationKit._getvectortype(::GridapProblem{Tfe, Tu}) where {Tfe, Tu} = Tu
BifurcationKit._getvectortype(pb::GridapBifProblem) = BifurcationKit._getvectortype(pb.probFE)
BifurcationKit.isinplace(::GridapBifProblem) = false
BifurcationKit.residual(pb::GridapBifProblem, u, p) = residual(pb.probFE, u, p)
BifurcationKit.jacobian(pb::GridapBifProblem, u, p) = jacobian(pb.probFE, u, p)
BifurcationKit.dF(pb::GridapBifProblem, u, p, dx) = BifurcationKit.apply(BifurcationKit.jacobian(pb, u, p), dx)
BifurcationKit.d2F(pb::GridapBifProblem, u, p, dx1, dx2) = pb.probFE(u, p, dx1, dx2)
BifurcationKit.d3F(pb::GridapBifProblem, u, p, dx1, dx2, dx3) = pb.probFE(u, p, dx1, dx2, dx3)
BifurcationKit.is_symmetric(::GridapBifProblem) = false
BifurcationKit.has_adjoint(::GridapBifProblem) = false
BifurcationKit.getdelta(pb::GridapBifProblem) = pb.δ
BifurcationKit.save_solution(::GridapBifProblem, x, p) = x
BifurcationKit.has_adjoint_MF(::GridapBifProblem) = false # TODO improve this using AD
BifurcationKit.update!(::GridapBifProblem, args...) = true
BifurcationKit.residual!(prob::GridapBifProblem, out, x, p) = out .= BK.residual(prob, x, p)

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
    probFE = GridapProblem(res, jacFE, d2res, d3res, V, U, nothing, dΩ, mass)
    # type unstable but simplifies the types a lot
    jet = BK.Jet(; δ = delta, R01, R02, R11, kwargs_jet...)
    return GridapBifProblem(probFE, Gridap.get_free_dof_values(u0), parms, lens, plot_solution, record_from_solution, delta, jet)
end

get_mass_matrix(prob::GridapBifProblem) = get_mass_matrix(prob.probFE)
BK.has_hessian(prob::GridapBifProblem) = BK.has_hessian(prob.VF)

BK.R01(prob::GridapBifProblem, x, p) = BK.R01(BK.has_R01_trait(prob.jet), prob, x, p)
BK.R01(::BK.TraitUserPassed, prob::GridapBifProblem, x, p) = prob.jet.R01(x, p)
BK.R02(prob::GridapBifProblem, x, p) = BK.R02(BK.has_R02_trait(prob.jet), prob, x, p)
BK.R02(::BK.TraitUserPassed, prob::GridapBifProblem, x, p) = prob.jet.R02(x, p)

BK.R11(prob::GridapBifProblem, x, p, dx) = BK.R11(BK.has_R11_trait(prob.jet), prob, x, p, dx)
BK.R11(::BK.TraitUserPassed, prob::GridapBifProblem, x, p, dx) = prob.jet.R11(x, p, dx)

BK.is_mass_matrix_constant(::GridapBifProblem) = true
BK.getmassmatrix(prob::GridapBifProblem, x, p) = get_mass_matrix(prob)