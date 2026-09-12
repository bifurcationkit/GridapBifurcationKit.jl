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

        @testset "mass derivatives (minimally augmented Hopf)" begin
            m = length(x1)
            v = collect(range(1.0, 2.0, length = m)); v ./= norm(v)
            w = collect(range(-0.5, 0.5, length = m)); w ./= norm(w)

            # ∇_x ⟨w, M(x, p) v⟩ by AutoDiff vs central finite differences
            g_ad = BifurcationKit.∇_x_mass_matrix(prob_sd, x1, par, v, w)
            g_fd = similar(g_ad)
            ϵ = 1e-6
            for i in eachindex(x1)
                ei = zero(x1); ei[i] = ϵ
                Mp = BifurcationKit.getmassmatrix(prob_sd, x1 + ei, par)
                Mm = BifurcationKit.getmassmatrix(prob_sd, x1 - ei, par)
                g_fd[i] = (BifurcationKit.dot_with_mass(w, Mp, v) -
                           BifurcationKit.dot_with_mass(w, Mm, v)) / (2ϵ)
            end
            @test g_ad ≈ g_fd rtol = 1e-5

            # ∂_p ⟨w, M(x, p) v⟩: λ-dependent SDMass, AutoDiff vs finite differences
            mass_sdp = (u, p, du, v) -> ∫(p.λ * (1 + u * u) * du * v) * dΩ
            prob_sdp = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                                        jac = jac, mass = mass_sdp)
            d_ad = BifurcationKit.R01_mass_matrix(prob_sdp, x1, par, v, w)
            ϵp = 1e-6
            Mp = BifurcationKit.getmassmatrix(prob_sdp, x1, (λ = par.λ + ϵp,))
            Mm = BifurcationKit.getmassmatrix(prob_sdp, x1, (λ = par.λ - ϵp,))
            d_fd = (BifurcationKit.dot_with_mass(w, Mp, v) -
                    BifurcationKit.dot_with_mass(w, Mm, v)) / (2ϵp)
            @test d_ad ≈ d_fd rtol = 1e-5

            # analytic state-derivative form `dM(u, p, du1, du2, v)` used for `∇xM`
            dM_sd = (u, p, du1, du2, v) -> ∫(2 * u * du1 * du2 * v) * dΩ
            prob_dM = GridapBifProblem(res, uh, par, V, U, dΩ, (@optic _.λ);
                                       jac = jac, mass = mass_sd, ∇xM = dM_sd)
            @test BifurcationKit.∇_x_mass_matrix(prob_dM, x1, par, v, w) ≈ g_ad rtol = 1e-9

            # complex directions (real / imaginary split)
            vc = v .+ im .* reverse(v)
            wc = w .+ im .* reverse(w)
            @test BifurcationKit.∇_x_mass_matrix(prob_dM, x1, par, vc, wc) ≈
                  BifurcationKit.∇_x_mass_matrix(prob_sd, x1, par, vc, wc) rtol = 1e-9
        end
    end

    # the AutoDiff mass gradient must also work when Dirichlet dofs are present
    # (mixing dual free dofs with Dirichlet values is rejected by Gridap).
    @testset "state-dependent mass with Dirichlet dofs" begin
        labels = get_face_labeling(model)
        add_tag_from_tags!(labels, "diri", [1, 2])
        Vd = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = ["diri"])
        Ud = TrialFESpace(Vd, 0.0)

        mass_sd = (u, p, du, v) -> ∫((1 + u * u) * du * v) * dΩ
        prob_d = GridapBifProblem(res, zero(Ud), par, Vd, Ud, dΩ, (@optic _.λ);
                                  jac = jac, mass = mass_sd)

        md = length(get_free_dof_values(zero(Ud)))
        xd = collect(range(0.1, 0.3, length = md))
        v = collect(range(1.0, 2.0, length = md)); v ./= norm(v)
        w = collect(range(-0.5, 0.5, length = md)); w ./= norm(w)

        g_ad = BifurcationKit.∇_x_mass_matrix(prob_d, xd, par, v, w)
        g_fd = similar(g_ad)
        ϵ = 1e-6
        for i in eachindex(xd)
            ei = zero(xd); ei[i] = ϵ
            Mp = BifurcationKit.getmassmatrix(prob_d, xd + ei, par)
            Mm = BifurcationKit.getmassmatrix(prob_d, xd - ei, par)
            g_fd[i] = (BifurcationKit.dot_with_mass(w, Mp, v) -
                       BifurcationKit.dot_with_mass(w, Mm, v)) / (2ϵ)
        end
        @test g_ad ≈ g_fd rtol = 1e-5
    end
end
