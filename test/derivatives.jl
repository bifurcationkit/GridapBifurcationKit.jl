using Test
using Gridap
using Gridap.FESpaces
using GridapBifurcationKit
using BifurcationKit

# Nonlinear 1D problem: -u'' + u + u³ = λ, so that the second and third
# derivatives of the residual are non trivial:
#   res(u, p, v)             = ∫(∇u⋅∇v + u v + u³ v - λ v)dΩ
#   jac(u, p, du, v)         = ∫(∇du⋅∇v + du v + 3 u² du v)dΩ
#   d2res(u, p, du1, du2, v)     = ∫(6 u du1 du2 v)dΩ
#   d3res(u, p, du1, du2, du3, v) = ∫(6 du1 du2 du3 v)dΩ
@testset "d2F / d3F derivative flavours" begin
    n = 8
    model = CartesianDiscreteModel((0.0, 1.0), (n,))
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = TestFESpace(model, reffe; conformity = :H1)
    U = TrialFESpace(V)
    dΩ = Measure(Triangulation(model), 2)

    res(u, p, v) = ∫(∇(u) ⋅ ∇(v) + u * v + u * u * u * v - p.λ * v) * dΩ
    jac(u, p, du, v) = ∫(∇(du) ⋅ ∇(v) + du * v + 3 * u * u * du * v) * dΩ
    d2res(u, p, du1, du2, v) = ∫(6 * u * du1 * du2 * v) * dΩ
    d3res(u, p, du1, du2, du3, v) = ∫(6 * du1 * du2 * du3 * v) * dΩ

    uh = zero(U)
    par = (λ = 1.0,)

    # the three flavours for `d2res`/`d3res`
    prob_an = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                               jac = jac, d2res = d2res, d3res = d3res)
    prob_ad = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                               jac = jac,
                               d2res = BifurcationKit.AutoDiff(),
                               d3res = BifurcationKit.AutoDiff())
    prob_fd = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                               jac = jac,
                               d2res = BifurcationKit.FiniteDifferences(),
                               d3res = BifurcationKit.FiniteDifferences())
    prob_nothing = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ); jac = jac)

    m = length(get_free_dof_values(uh))
    x = collect(range(-0.3, 0.4, length = m))
    du1 = collect(range(0.5, 1.0, length = m))
    du2 = collect(range(-1.0, -0.5, length = m))
    du3 = collect(range(0.0, 1.0, length = m))

    v2_an = prob_an.probFE(x, par, du1, du2)
    v2_ad = prob_ad.probFE(x, par, du1, du2)
    v2_fd = prob_fd.probFE(x, par, du1, du2)
    v2_no = prob_nothing.probFE(x, par, du1, du2)

    v3_an = prob_an.probFE(x, par, du1, du2, du3)
    v3_ad = prob_ad.probFE(x, par, du1, du2, du3)
    v3_fd = prob_fd.probFE(x, par, du1, du2, du3)

    # sanity: the derivatives are non trivial
    @test norm(v2_an) > 0
    @test norm(v3_an) > 0

    # AutoDiff and FiniteDifferences reproduce the analytic forms
    @test v2_ad ≈ v2_an rtol = 1e-9
    @test v3_ad ≈ v3_an rtol = 1e-9
    @test v2_fd ≈ v2_an rtol = 1e-5
    @test v3_fd ≈ v3_an rtol = 1e-4

    # `nothing` is the finite-difference default
    @test v2_no ≈ v2_fd rtol = 1e-8

    # the BifurcationKit `d2F` / `d3F` wrappers agree as well
    @test BifurcationKit.d2F(prob_ad, x, par, du1, du2) ≈ v2_ad rtol = 1e-10
    @test BifurcationKit.d3F(prob_ad, x, par, du1, du2, du3) ≈ v3_ad rtol = 1e-10
end

# The `AutoDiff` path must also work for a multi-field problem with Dirichlet
# boundary conditions (the situation of the incompressible flow examples).
@testset "d2F / d3F AutoDiff on a multi-field problem" begin
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (4, 4))
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "diri", collect(1:8))

    reffe_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
    reffe_p = ReferenceFE(lagrangian, Float64, 1)
    V = TestFESpace(model, reffe_u; conformity = :H1, dirichlet_tags = ["diri"])
    U = TrialFESpace(V, VectorValue(0.0, 0.0))
    Q = TestFESpace(model, reffe_p; conformity = :L2)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U, P])
    Y = MultiFieldFESpace([V, Q])
    dΩ = Measure(Triangulation(model), 2)

    res((u, p), par, (v, q)) = ∫(∇(u) ⊙ ∇(v) + u ⋅ v + (u ⋅ u) * (u ⋅ v) +
                                 p * (∇ ⋅ v) + q * (∇ ⋅ u) -
                                 par.λ * (v ⋅ VectorValue(1.0, 0.0))) * dΩ
    jac((u, p), par, (du, dp), (v, q)) = ∫(∇(du) ⊙ ∇(v) + du ⋅ v +
                                           2 * (u ⋅ du) * (u ⋅ v) + (u ⋅ u) * (du ⋅ v) +
                                           dp * (∇ ⋅ v) + q * (∇ ⋅ du)) * dΩ

    uh = zero(X)
    par = (λ = 1.0,)
    prob_ad = GridapBifProblem(res, uh, par, Y, X, dΩ, (@optic _.λ);
                               jac = jac,
                               d2res = BifurcationKit.AutoDiff(),
                               d3res = BifurcationKit.AutoDiff())
    prob_fd = GridapBifProblem(res, uh, par, Y, X, dΩ, (@optic _.λ);
                               jac = jac,
                               d2res = BifurcationKit.FiniteDifferences(),
                               d3res = BifurcationKit.FiniteDifferences())

    m = length(get_free_dof_values(uh))
    x = collect(range(-0.2, 0.3, length = m))
    du1 = collect(range(-0.4, 0.5, length = m))
    du2 = collect(range(0.6, 1.1, length = m))
    du3 = collect(range(-1.0, -0.2, length = m))

    v2_ad = prob_ad.probFE(x, par, du1, du2)
    v2_fd = prob_fd.probFE(x, par, du1, du2)
    v3_ad = prob_ad.probFE(x, par, du1, du2, du3)
    v3_fd = prob_fd.probFE(x, par, du1, du2, du3)

    @test norm(v2_ad) > 0
    @test v2_ad ≈ v2_fd rtol = 1e-6
    @test v3_ad ≈ v3_fd rtol = 1e-4
end
