# example/cylinder_defs_kernel.jl
#
# Drop-in replacement for cylinder_defs.jl in which `res` and `jac` are computed
# by a hand-written cell kernel instead of a tree of CellField operations.
# Shared definitions for the 2-D compressible flow past a cylinder example
# (Section 3.3 of Douglas & Jolivet, CPC 326 (2026) 110221).
#
# Unknowns q = (u, T, p): u is the 2-D velocity (P2), T the temperature (P2) and
# p the pressure (P1).  Parameters are (1/Re, 1/Pr, Ma^2, gamma).  The mass
# operator M(q) is state-dependent (a DAE `M(q) q' = F(q)`).

using Gridap

_u_zero(x) = VectorValue(0.0, 0.0)
_u_inlet(x) = VectorValue(1.0, 0.0)

cyl_par(; Re = 40.0, M = 0.0, Pr = 0.72, γ = 1.4) =
    (;Re, Pr, M, γ)

"""
    cylinder_spaces(model; full=false)

Mixed Taylor-Hood spaces `(X, Y) = (trial, test)` on `model`. With `full=false`
(the symmetric half) the tag `"sym"` carries `u_y = 0`; with `full=true` the
mirrored full domain has no symmetry axis.  Boundary tags: `"wall"`, `"inlet"`,
`"lat"` (+ `"sym"` for the half).
"""
function cylinder_spaces(model; full = false)
    reffe_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, 2)
    reffe_T = ReferenceFE(lagrangian, Float64, 2)
    reffe_p = ReferenceFE(lagrangian, Float64, 1)

    tags = full ? ["wall", "inlet", "lat"] : ["wall", "inlet", "lat", "sym"]
    masks = full ? [(true,true), (true,true), (false,true)] :
                   [(true,true), (true,true), (false,true), (false,true)]

    V_u = TestFESpace(model, reffe_u; conformity = :H1,
                      dirichlet_tags = tags, dirichlet_masks = masks)
    V_T = TestFESpace(model, reffe_T; conformity = :H1, dirichlet_tags = ["inlet"])
    V_p = TestFESpace(model, reffe_p; conformity = :H1)

    uvals = Any[t == "inlet" ? (_u_inlet) : (_u_zero) for t in tags]
    U_u = TrialFESpace(V_u, uvals)
    U_T = TrialFESpace(V_T, x -> 1.0)
    U_p = TrialFESpace(V_p)

    X = MultiFieldFESpace([U_u, U_T, U_p])
    Y = MultiFieldFESpace([V_u, V_T, V_p])
    return X, Y
end

# distance to physical domain
bsds_fun(x) = 1.0e-4 * (max(x[1] - 50, -20 - x[1], 0)^2 +
                        max(abs(x[2]) - 20, 0)^2)

"""
    cylinder_measures(model)

Return `(Ω, dΩ, bsds)` where `bsds` is the absorbing-layer weight
`1e-4 (max(x-50, -20-x, 0)^2 + max(|y|-20, 0)^2)` as a `CellField`.
"""
function cylinder_measures(model, degree)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    bsds = CellField(bsds_fun, Ω)
    return Ω, dΩ, bsds
end

ugradf(U, f) = U ⋅ ∇(f)
diffg(g, f) = ∇(g) ⋅ ∇(f)
vdotu(v, 𝐮) = v ⋅ 𝐮
divu(𝐮) = tr(∇(𝐮))
ugradu(v, U, 𝐮) = v ⋅ (∇(𝐮)' ⋅ U)

τ(𝐮) = ∇(𝐮) + ∇(𝐮)' - (2 / 3) * divu(𝐮) * one(TensorValue{2,2,Float64})
visc(v, 𝐮) = ∇(v) ⊙ τ(𝐮)
visc2(v, T, 𝐮) = (v ⊗ ∇(T)) ⊙ τ(𝐮)

using Gridap.Geometry, Gridap.Arrays, Gridap.Fields, Gridap.CellData
using Gridap.Fields: meas
import Gridap.Arrays: return_cache, evaluate!

# Pointwise physics.  The residual is written per test function,
#   r = ∫ v⋅Fu + ∇v⊙Gu + vT*FT + ∇vT⋅GT + vp*Fp,
# and the Jacobian as the directional derivatives of (Fu, Gu, FT, GT, Fp).
const _I2 = one(TensorValue{2,2,Float64})
const _ex = VectorValue(1.0, 0.0)

function kernel_coeffs(pars)
    invRe = 1 / pars.Re; invPr = 1 / pars.Pr; γ = pars.γ; Ma2 = pars.M^2
    a = γ * Ma2
    (; a, γ, Ma2, invRe, ediff = γ * invRe * invPr, eviscl = a * (γ - 1) * invRe, σ = -1.0)
end

@inline τ_p(∇u) = ∇u + ∇u' - (2 / 3) * tr(∇u) * _I2

@inline function res_flux(c, u, ∇u, T, ∇T, p, ∇p, s)
    ap = c.a * p + 1
    τu = τ_p(∇u)
    Fu = ap * (∇u' ⋅ u) - p * ∇T + c.invRe * (τu ⋅ ∇T) + s * (u - _ex)
    Gu = (c.invRe * T) * τu - (T * p) * _I2
    FT = ap * (u ⋅ ∇T + (c.γ - 1) * T * tr(∇u)) - c.eviscl * T * (∇u ⊙ τu) +
         c.ediff * (∇T ⋅ ∇T) + (s * c.γ) * (T - 1)
    GT = (c.ediff * T) * ∇T
    Fp = c.a * T * (u ⋅ ∇p) - ap * (u ⋅ ∇T - T * tr(∇u)) + (s * c.γ * c.Ma2) * p
    return Fu, Gu, FT, GT, Fp
end

@inline function jac_flux_u(c, u, ∇u, T, ∇T, p, ∇p, s, δu, ∇δu)
    ap = c.a * p + 1
    τu = τ_p(∇u); τδ = τ_p(∇δu)
    dFu = ap * (∇δu' ⋅ u + ∇u' ⋅ δu) + c.invRe * (τδ ⋅ ∇T) + s * δu
    dGu = (c.invRe * T) * τδ
    dFT = ap * (δu ⋅ ∇T + (c.γ - 1) * T * tr(∇δu)) - c.eviscl * T * (∇δu ⊙ τu + ∇u ⊙ τδ)
    dFp = c.a * T * (δu ⋅ ∇p) - ap * (δu ⋅ ∇T - T * tr(∇δu))
    return dFu, dGu, dFT, dFp
end

@inline function jac_flux_T(c, u, ∇u, T, ∇T, p, ∇p, s, δT, ∇δT)
    ap = c.a * p + 1
    τu = τ_p(∇u)
    dFu = -p * ∇δT + c.invRe * (τu ⋅ ∇δT)
    dGu = (c.invRe * δT) * τu - (δT * p) * _I2
    dFT = ap * (u ⋅ ∇δT + (c.γ - 1) * δT * tr(∇u)) - c.eviscl * δT * (∇u ⊙ τu) +
          2 * c.ediff * (∇T ⋅ ∇δT) + (s * c.γ) * δT
    dGT = c.ediff * (δT * ∇T + T * ∇δT)
    dFp = c.a * δT * (u ⋅ ∇p) - ap * (u ⋅ ∇δT - δT * tr(∇u))
    return dFu, dGu, dFT, dGT, dFp
end

@inline function jac_flux_p(c, u, ∇u, T, ∇T, p, ∇p, s, δp, ∇δp)
    dFu = (c.a * δp) * (∇u' ⋅ u) - δp * ∇T
    dGu = -(T * δp) * _I2
    dFT = (c.a * δp) * (u ⋅ ∇T + (c.γ - 1) * T * tr(∇u))
    dFp = c.a * T * (u ⋅ ∇δp) - (c.a * δp) * (u ⋅ ∇T - T * tr(∇u)) + (s * c.γ * c.Ma2) * δp
    return dFu, dGu, dFT, dFp
end

# Cell kernels.  Arguments, per cell: quadrature weights w, cell-map Jacobian
# Jt at the points, basis values and gradients (np × ndof) of the three
# fields, then the state and the sponge weight at the points.
struct ResKernel{C} <: Map
    c::C
end
struct JacKernel{C} <: Map
    c::C
end

_ndofs(φu, φT, φp) = (size(φu, 2), size(φT, 2), size(φp, 2))

# scalar entry type of a (possibly vector/tensor-valued) cell data type, e.g.
# `VectorValue{Dual} -> Dual`; used to type the kernel caches with the (possibly
# dual) eltype of the state, as required by the AD-based derivatives.
_scalar_eltype(::Type{T}) where {T} = (U = eltype(T); U === T ? T : _scalar_eltype(U))

# eltype of the kernel caches: promoted over the state fields. With Gridap's
# split AD only one field is dual at a time, but all residual blocks can
# receive dual contributions through the cross-coupling terms.
_kernel_eltype(u, T, p) = promote_type(_scalar_eltype(eltype(u)), eltype(T), eltype(p))

function return_cache(::ResKernel, w, Jt, φu, ∇φu, φT, ∇φT, φp, ∇φp, u, ∇u, T, ∇T, p, ∇p, s)
    nu, nT, np = _ndofs(φu, φT, φp)
    𝒯 = _kernel_eltype(u, T, p)
    ArrayBlock([zeros(𝒯, nu), zeros(𝒯, nT), zeros(𝒯, np)], [true, true, true])
end

function evaluate!(r, k::ResKernel, w, Jt, φu, ∇φu, φT, ∇φT, φp, ∇φp, u, ∇u, T, ∇T, p, ∇p, s)
    c = k.c
    nu, nT, np = _ndofs(φu, φT, φp)
    ru, rT, rp = r.array
    (length(ru), length(rT), length(rp)) == (nu, nT, np) || error("non-uniform cell dofs")
    fill!(ru, 0.0); fill!(rT, 0.0); fill!(rp, 0.0)
    for q in eachindex(w)
        dV = c.σ * w[q] * meas(Jt[q])
        Fu, Gu, FT, GT, Fp = res_flux(c, u[q], ∇u[q], T[q], ∇T[q], p[q], ∇p[q], s[q])
        for i in 1:nu
            ru[i] += dV * (φu[q, i] ⋅ Fu + ∇φu[q, i] ⊙ Gu)
        end
        for i in 1:nT
            rT[i] += dV * (φT[q, i] * FT + ∇φT[q, i] ⋅ GT)
        end
        for i in 1:np
            rp[i] += dV * (φp[q, i] * Fp)
        end
    end
    return r
end

function return_cache(::JacKernel, w, Jt, φu, ∇φu, φT, ∇φT, φp, ∇φp, u, ∇u, T, ∇T, p, ∇p, s)
    n = _ndofs(φu, φT, φp)
    𝒯 = _kernel_eltype(u, T, p)
    ArrayBlock([zeros(𝒯, n[i], n[j]) for i in 1:3, j in 1:3], fill(true, 3, 3))
end

function evaluate!(J, k::JacKernel, w, Jt, φu, ∇φu, φT, ∇φT, φp, ∇φp, u, ∇u, T, ∇T, p, ∇p, s)
    c = k.c
    nu, nT, np = _ndofs(φu, φT, φp)
    size(J.array[1, 1]) == (nu, nu) && size(J.array[2, 2]) == (nT, nT) && size(J.array[3, 3]) == (np, np) ||
        error("non-uniform cell dofs")
    A = J.array
    foreach(a -> fill!(a, 0.0), A)
    for q in eachindex(w)
        dV = c.σ * w[q] * meas(Jt[q])
        st = (u[q], ∇u[q], T[q], ∇T[q], p[q], ∇p[q], s[q])
        for j in 1:nu
            dFu, dGu, dFT, dFp = jac_flux_u(c, st..., φu[q, j], ∇φu[q, j])
            for i in 1:nu; A[1, 1][i, j] += dV * (φu[q, i] ⋅ dFu + ∇φu[q, i] ⊙ dGu); end
            for i in 1:nT; A[2, 1][i, j] += dV * (φT[q, i] * dFT); end
            for i in 1:np; A[3, 1][i, j] += dV * (φp[q, i] * dFp); end
        end
        for j in 1:nT
            dFu, dGu, dFT, dGT, dFp = jac_flux_T(c, st..., φT[q, j], ∇φT[q, j])
            for i in 1:nu; A[1, 2][i, j] += dV * (φu[q, i] ⋅ dFu + ∇φu[q, i] ⊙ dGu); end
            for i in 1:nT; A[2, 2][i, j] += dV * (φT[q, i] * dFT + ∇φT[q, i] ⋅ dGT); end
            for i in 1:np; A[3, 2][i, j] += dV * (φp[q, i] * dFp); end
        end
        for j in 1:np
            dFu, dGu, dFT, dFp = jac_flux_p(c, st..., φp[q, j], ∇φp[q, j])
            for i in 1:nu; A[1, 3][i, j] += dV * (φu[q, i] ⋅ dFu + ∇φu[q, i] ⊙ dGu); end
            for i in 1:nT; A[2, 3][i, j] += dV * (φT[q, i] * dFT); end
            for i in 1:np; A[3, 3][i, j] += dV * (φp[q, i] * dFp); end
        end
    end
    return J
end

"""
    cylinder_forms(bsds, dΩ)

Return `(res, jac, mass, dM, d2res, d3res, R01_Re, R11_Re)` translated from
`eqns_douglas_jolivet_2026_cylinder.idp`. `res(u,p,v)`, `jac(u,p,du,v)`, the
state-dependent `mass(u,p,du,v)`, its state derivative
`dM(u,p,du1,du2,v)` (the `vdM` operator, to be passed as `∇xM`), the
Hessian/Tressian of the residual `d2res(u,p,du1,du2,v)`,
`d3res(u,p,du1,du2,du3,v)` (the latter translated from the stationary,
`iomega = 0`, parts of `vH`/`vT`), and the parameter derivatives w.r.t. `Re`
`R01_Re(u,p,v)` (∂_Re res) and `R11_Re(u,p,du,v)` (∂²res/∂Re∂u = ∂_Re jac).
"""
function cylinder_forms(bsds, dΩ)
    quad = dΩ.quad
    Ω = quad.trian
    xq = get_cell_points(quad)
    w = quad.cell_weight
    Jt = lazy_map(evaluate, lazy_map(∇, get_cell_map(Ω)), quad.cell_point)
    sq = bsds(xq)
    # Trial and test shape functions coincide; the test basis serves for both.
    basis((v, vT, vp)) = (v.single_field(xq), ∇(v.single_field)(xq),
                          vT.single_field(xq), ∇(vT.single_field)(xq),
                          vp.single_field(xq), ∇(vp.single_field)(xq))
    state((𝐮, 𝐓, 𝐩)) = (𝐮(xq), ∇(𝐮)(xq), 𝐓(xq), ∇(𝐓)(xq), 𝐩(xq), ∇(𝐩)(xq), sq)
    contribution(cell) = add_contribution!(DomainContribution(), Ω, cell)

    res = function (x, pars, y)
        contribution(lazy_map(ResKernel(kernel_coeffs(pars)), w, Jt, basis(y)..., state(x)...))
    end

    jac = function (x, pars, dx, y)
        contribution(lazy_map(JacKernel(kernel_coeffs(pars)), w, Jt, basis(y)..., state(x)...))
    end

    R01_Re = function ((u, T, pr), pars, (v, vT, vp))
        @error "R01_res"
        σ = -1.0
        invRe  = 1 / pars.Re
        invRe2 = -invRe^2                      # ∂_{Re}(1/Re)
        a      = pars.γ * pars.M^2
        invPr  = 1 / pars.Pr
        return σ * (∫(
            invRe2 * ( (T * visc(v, u) + visc2(v, T, u))
                    - vT * a * (pars.γ - 1) * T * visc(u, u)
                    + pars.γ * invPr * (T * diffg(vT, T) + vT * diffg(T, T)) )
        ) * dΩ)
    end

    # mixed derivative ∂²res/∂Re∂u, i.e. ∂_{Re} jac: a *bilinear* form in
    # (du, v) (used by `BK.R11`). Like `R01_Re`, only the `1/Re` factors of
    # `jac` contribute, so the bracket is exactly the `invRe`-multiplied part
    # of `jac`.
    R11_Re = function ((u, T, pr), pars, (du, duT, dup), (v, vT, vp))
        @error "R11_Re"
        σ = -1.0
        invRe  = 1 / pars.Re
        invRe2 = -invRe^2                      # ∂_{Re}(1/Re)
        a      = pars.γ * pars.M^2
        invPr  = 1 / pars.Pr
        return σ * (∫(
            invRe2 * ( duT * visc(v, u) + T * visc(v, du)
                     + visc2(v, duT, u) + visc2(v, T, du)
                     - vT * a * (pars.γ - 1) *
                           (duT * visc(u, u) + T * (visc(du, u) + visc(u, du)))
                     + pars.γ * invPr * (duT * diffg(vT, T) + T * diffg(vT, duT)
                                         + 2.0 * vT * diffg(T, duT)) )
        ) * dΩ)
    end

    mass = function ((u, T, pr), pars, (du, duT, dup), (v, vT, vp))
        a = pars.γ * pars.M^2
        return ∫((a * pr + 1) * (vdotu(v, du) + (vT - vp) * duT)
                 + a * T * vp * dup) * dΩ
    end

    # matrix-free application M(x, p) * dx: weak form of the mass bilinear
    # operator (state `u`, direction `du`), passed as `applyM` to
    # `GridapBifProblem`, which assembles it into a free-dof vector.
    apply_mass = function ((u, T, pr), pars, (du, duT, dup), (v, vT, vp))
        a = pars.γ * pars.M^2
        return ∫((a * pr + 1) * (vdotu(v, du) + (vT - vp) * duT)
                 + a * T * vp * dup) * dΩ
    end

    # MASS DERIVATIVE OPERATOR
    # Trilinear form ∂_u⟨(v, vT, vp), M(u, p)(du2, du2T, du2p)⟩ in the state
    # direction (du1, du1T, du1p). Translation of `vdM` from
    # plays the role of the mass direction of `∇_x_mass_matrix`.
    dM = function ((u, T, pr), pars, (du1, du1T, du1p), (du2, du2T, du2p), (v, vT, vp))
        a = pars.γ * pars.M^2
        return ∫(a * (du1p * (vdotu(v, du2) + (vT - vp) * du2T)
                      + du1T * vp * du2p)) * dΩ
    end

    # Hessian (second derivative of the residual), translated from the
    # iomega = 0 part of the first block of `vH` in
    d2res = function ((u, T, pr), pars, (du1, du1T, du1p), (du2, du2T, du2p), (v, vT, vp))
        @error "d2res"
        invRe = 1 / pars.Re
        invPr = 1 / pars.Pr
        a = pars.γ * pars.M^2
        σ = -1.0
        return σ * (∫(
            a * (du1p * (ugradu(v, du2, u) + ugradu(v, u, du2))
                 + du2p * (ugradu(v, du1, u) + ugradu(v, u, du1)))
            + (a * pr + 1) * (ugradu(v, du1, du2) + ugradu(v, du2, du1))
            - (du2T * divu(v) + ugradf(v, du2T)) * du1p
            - (du1T * divu(v) + ugradf(v, du1T)) * du2p
            + invRe * (du1T * visc(v, du2) + du2T * visc(v, du1)
                       + visc2(v, du1T, du2) + visc2(v, du2T, du1))
            + vT * (
                a * (du1p * (ugradf(du2, T) + ugradf(u, du2T) + (pars.γ - 1) * (du2T * divu(u) + T * divu(du2)))
                     + du2p * (ugradf(du1, T) + ugradf(u, du1T) + (pars.γ - 1) * (du1T * divu(u) + T * divu(du1))))
                + (a * pr + 1) * (ugradf(du1, du2T) + ugradf(du2, du1T)
                                  + (pars.γ - 1) * (du1T * divu(du2) + du2T * divu(du1)))
                - a * ((pars.γ - 1) * invRe * (du1T * (visc(du2, u) + visc(u, du2))
                                                + du2T * (visc(du1, u) + visc(u, du1))
                                                + T * (visc(du1, du2) + visc(du2, du1))))
            )
            + pars.γ * invRe * invPr * (du1T * diffg(vT, du2T) + du2T * diffg(vT, du1T)
                                            + 2.0 * vT * diffg(du2T, du1T))
            + vp * (
                  a * (du1T * (ugradf(du2, pr) + ugradf(u, du2p))
                     + du2T * (ugradf(du1, pr) + ugradf(u, du1p))
                     + T * (ugradf(du1, du2p) + ugradf(du2, du1p))
                     + du1p * (du2T * divu(u) + T * divu(du2) - ugradf(du2, T) - ugradf(u, du2T))
                     + du2p * (du1T * divu(u) + T * divu(du1) - ugradf(du1, T) - ugradf(u, du1T)))
                + (a * pr + 1) * (du1T * divu(du2) + du2T * divu(du1)
                                  - ugradf(du1, du2T) - ugradf(du2, du1T))
            )
        ) * dΩ)
    end

    # Tressian (third derivative of the residual), translated from the first
    # It is constant in the state (the residual is cubic).
    d3res = function ((u, T, pr), pars, (du1, du1T, du1p), (du2, du2T, du2p), (du3, du3T, du3p), (v, vT, vp))
        invRe = 1 / pars.Re
        a = pars.γ * pars.M^2
        σ = -1.0
        return σ * (∫(
              a * (du1p * (ugradu(v, du2, du3) + ugradu(v, du3, du2))
                 + du2p * (ugradu(v, du1, du3) + ugradu(v, du3, du1))
                 + du3p * (ugradu(v, du1, du2) + ugradu(v, du2, du1)))
            + vT * (
                  a * (du1p * (ugradf(du2, du3T) + ugradf(du3, du2T) + (pars.γ - 1) * (du2T * divu(du3) + du3T * divu(du2)))
                     + du2p * (ugradf(du1, du3T) + ugradf(du3, du1T) + (pars.γ - 1) * (du1T * divu(du3) + du3T * divu(du1)))
                     + du3p * (ugradf(du1, du2T) + ugradf(du2, du1T) + (pars.γ - 1) * (du1T * divu(du2) + du2T * divu(du1))))
                   - a * ((pars.γ - 1) * invRe * (du1T * (visc(du2, du3) + visc(du3, du2))
                                                + du2T * (visc(du1, du3) + visc(du3, du1))
                                                + du3T * (visc(du1, du2) + visc(du2, du1))))
            )
            + vp * a * (
                  du1T * (ugradf(du2, du3p) + ugradf(du3, du2p))
                + du2T * (ugradf(du1, du3p) + ugradf(du3, du1p))
                + du3T * (ugradf(du1, du2p) + ugradf(du2, du1p))
                + du1p * (du2T * divu(du3) + du3T * divu(du2) - ugradf(du2, du3T) - ugradf(du3, du2T))
                + du2p * (du1T * divu(du3) + du3T * divu(du1) - ugradf(du1, du3T) - ugradf(du3, du1T))
                + du3p * (du1T * divu(du2) + du2T * divu(du1) - ugradf(du1, du2T) - ugradf(du2, du1T))
            )
        ) * dΩ)
    end

    return (;res, jac, mass, dM, d2res, d3res, R01_Re, R11_Re, apply_mass)
end

"""
    prep_cylinder_drag(model, order; tags = ["wall"], degree = 2*order)

Prepare the cylinder `wall` boundary (boundary triangulation, measure, outward
normal) and return the closure `drag(q, pars) -> (; cd, Fx)`.

- `q` is the current state (`MultiFieldFEFunction` with `q[1]=u`, `q[2]=T`,
  `q[3]=pr`), e.g. `qrec` after `qrec.free_values .= x`;
- `pars` is the current value of the continuation parameter, i.e. `Re` (the
  integral uses `1/Re`).

Faithful translation of the `Fx` monitor of
`settings_douglas_jolivet_2026_cylinder.idp:48` (`1/Re` only, no `T`, viscous
traction `(∇u + ∇uᵀ − (∇·u) I)·n`), with the 2-D convention `cd = 2 Fx`
(`U∞ = D = ρ = 1`).
"""
function prep_cylinder_drag(model, order; tags = ["wall"], degree = 2 * order)
    ∂Ω = BoundaryTriangulation(model; tags = tags)
    dΓ = Measure(∂Ω, degree)
    n  = get_normal_vector(∂Ω)
    ex = VectorValue(1.0, 0.0)
    function drag(q, Re)
        u, T, pr = q[1], q[2], q[3]
        invRe = 1 / Re
        tv = (∇(u) + ∇(u)') ⋅ n - divu(u) * n
        Fx = sum(∫( pr * (n ⋅ ex) - invRe * (tv ⋅ ex) ) * dΓ)   # <- sum -> scalaire
        return (; cd = 2.0 * Fx, Fx)
    end
    return drag
end
