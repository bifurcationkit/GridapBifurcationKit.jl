using Test
using Gridap
using Gridap.FESpaces
using GridapBifurcationKit
using BifurcationKit

# simple 1D FE problem: -u'' + u = λ on (0,1)
# res(u, p, v) = ∫(∇u⋅∇v + u v)dΩ - λ∫v dΩ
# jac(u, p, du, v) = ∫(∇du⋅∇v + du v)dΩ
@testset "GridapBifProblem" begin
    n = 8
    model = CartesianDiscreteModel((0.0, 1.0), (n,))
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = TestFESpace(model, reffe, conformity = :H1)
    U = TrialFESpace(V)
    dΩ = Measure(Triangulation(model), 2)

    res(u, p, v) = ∫(∇(u) ⋅ ∇(v) + u * v - p.λ * v) * dΩ
    jac(u, p, du, v) = ∫(∇(du) ⋅ ∇(v) + du * v) * dΩ

    uh = zero(U)
    par = (λ = 1.0,)
    prob = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ); jac = jac)

    x = get_free_dof_values(uh)

    # zero state is in the kernel of the (linear) jacobian
    @test norm(BifurcationKit.jacobian(prob, x, par) * x) ≈ 0 atol = 1e-12

    # the source term makes the residual non trivial
    r = BifurcationKit.residual(prob, x, par)
    @test norm(r) > 0

    # the problem is self-adjoint: the jacobian is symmetric
    J = BifurcationKit.jacobian(prob, x, par)
    @test norm(J - transpose(J)) ≈ 0 atol = 1e-12

    # dF is the jacobian-vector product, consistent with finite differences of res
    dx = collect(range(0.0, 1.0, length = length(x)))
    ϵ = 1e-6
    fd = (BifurcationKit.residual(prob, x + ϵ * dx, par) - r) / ϵ
    @test fd ≈ BifurcationKit.dF(prob, x, par, dx) rtol = 1e-5

    @testset "mass interface" begin
        # 1. MassDefaut: no `mass` keyword, L² mass ∫(u⋅v)dΩ
        M = GridapBifurcationKit.get_mass_matrix(prob)
        m = length(x)
        @test size(M) == (m, m)
        @test GridapBifurcationKit.get_mass_matrix(prob, x, par) ≈ M
        @test BifurcationKit.getmassmatrix(prob, x, par) ≈ M
        @test BifurcationKit.is_mass_matrix_constant(prob)

        # 2. ConstMass: state-independent integrand mass(u, v)
        c = 3.0
        mass_c = (u, v) -> ∫(c * u * v) * dΩ
        prob_c = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                                  jac = jac, mass = mass_c)
        @test BifurcationKit.is_mass_matrix_constant(prob_c)
        @test GridapBifurcationKit.get_mass_matrix(prob_c) ≈ c * M
        @test GridapBifurcationKit.get_mass_matrix(prob_c, x, par) ≈ c * M
        @test BifurcationKit.getmassmatrix(prob_c, x, par) ≈ c * M

        # 3. SDMass: state-dependent integrand mass(u, p, du, v), i.e. M(x, p)
        mass_sd = (u, p, du, v) -> ∫((1 + u * u) * du * v) * dΩ
        prob_sd = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                                   jac = jac, mass = mass_sd)
        @test !BifurcationKit.is_mass_matrix_constant(prob_sd)
        # the state-less interface must refuse a state-dependent mass
        @test_throws ArgumentError GridapBifurcationKit.get_mass_matrix(prob_sd)
        # at the zero state the operator reduces to the default L² mass
        @test GridapBifurcationKit.get_mass_matrix(prob_sd, x, par) ≈ M
        # at a non trivial state it genuinely depends on the state
        u1 = interpolate(z -> z[1], U)
        x1 = get_free_dof_values(u1)
        M1 = GridapBifurcationKit.get_mass_matrix(prob_sd, x1, par)
        @test norm(M1 - M) > 0
        @test BifurcationKit.getmassmatrix(prob_sd, x1, par) ≈ M1
    end
end
