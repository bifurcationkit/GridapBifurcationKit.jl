using Revise, GLMakie
Makie.inline!(true)
using Gridap, SparseArrays, BifurcationKit
# using KrylovKit, LinearAlgebra
using Gridap.FESpaces
using GridapBifurcationKit
using GridapMakie
const BK = BifurcationKit
#############################################
# D = 6.5
# const par_seco = (D = D, τ = 0.05, α = 0.035, j0 = 1.2)
const par_seco = (D = 8., τ = 0.05, α = 0.02, j0 = 3.6)

begin
Nx = 201
const lx = 8pi / (par_seco.α*par_seco.τ/par_seco.D)^(1/4)
domain = (-lx, lx)
cells = (Nx)
model = CartesianDiscreteModel(domain, cells)

order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V0 = TestFESpace(model, reffe, conformity=:H1)
U0 = TrialFESpace(V0)

Y = MultiFieldFESpace([V0, V0])
const X = MultiFieldFESpace([U0, U0])

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω, degree)
end

function res((u1,u2),par,(V1, V2))
    (;D, τ, α, j0) = par
    NL1(u1,u2) =  (u2-u1)/((u2-u1)^2+1) - τ * u1
    NL2(u1,u2) = α * (j0 - (u2 - u1))

    ∫(  -∇(u1)⋅∇(V1) + NL1∘(u1,u2) ⋅ V1 +
    (-D) * ∇(u2)⋅∇(V2) + NL2∘(u1,u2) ⋅ V2 )*dΩ
end

function jac((u,v),par,(du, dv),(V1, V2))
    (;D, τ, α, j0) = par
    d1NL1(u1,u2) = -1/((u2 - u1)^2 + 1) - (u2 - u1)*(-2*u2 + 2*u1)/((u2 - u1)^2 + 1)^2 - τ
    d2NL1(u1,u2) = 1/((u2 - u1)^2 + 1) - (u2 - u1)*(2*u2 - 2*u1)/((u2 - u1)^2 + 1)^2
    d1NL2(u1,u2) = α
    d2NL2(u1,u2) = -α

    ∫(  -∇(du)⋅∇(V1) + (-D) * ∇(dv)⋅∇(V2) +
            d1NL1∘(u,v)⋅V1⋅du + d2NL1∘(u,v)⋅V1⋅dv +
            d1NL2∘(u,v)⋅V2⋅du + d2NL2∘(u,v)⋅V2⋅dv )*dΩ
end

uh = zero(X)
uh.free_values .= 1
####################################################################################################
function ustar(p)
    par = @set par_seco.j0 = p
    (;D, τ, α, j0) = par
    j0 /(j0^2+1)/τ
end
####################################################################################################
N2(u) = sqrt((∫(u⋅u)*dΩ |> sum) / (∫(1)*dΩ |> sum))

DOT(u,v) = (∫(u⋅v)*dΩ |> sum)
DOT(uh,uh) |> sum
function DOT(x::Vector, y::Vector)
    ux = zero(X)
    uy = zero(X)
    ux.free_values .= x
    uy.free_values .= y
    return DOT(ux, uy)/2/lx
end
####################################################################################################
function N2(x::AbstractVector, _X, us)
    uh = zero(X)
    uh.free_values .= x .- us
    N2(uh[1])
end

function NormL2T(prob_po, po::AbstractArray{Tu}, us) where {Tu}
    solt = BK.get_periodic_orbit(prob_po, po, nothing)
    u = solt.u

    nrm = zero(Tu)
    n = size(u, 2)

    for j in axes(u, 2)
        nrm += N2(u[:, j], X, us)^2
    end
    s = 1
    sqrt(nrm/(n * s))
end

recordSolSECO(x, p; k...) = (
            nrm =  N2(x, X, ustar(p)),
            norminf = norm(x[1:length(x)÷2] .- ustar(p), Inf),
            n2 = norm(x[1:length(x)÷2] .- ustar(p), 2))

plotSECO(ax,x,p;k...) = plotsol!(ax,x; k...)

function plotsol!(ax, sol; nd = 400, k...)
    solr, soli = sol
    lines!(ax, Ω, solr; label = "", k...)
end

function plotsol!(ax, x::AbstractVector; ax1 = nothing, iter = nothing, state = nothing, i = 1, k...)
    uh = zero(X)
    uh.free_values .= x
    plotsol!(ax, uh; k...)
end
plotsol(sol; k...) = (f=Figure();ax=Axis(f[1,1]);plotsol!(ax, sol; k...);f)

prob = GridapBifProblem(res, uh, par_seco, Y, X, dΩ, (@optic _.j0);
            jac = jac,
            jacobian_type = BK.FullSparseInplace(),
            plot_solution = plotSECO,
            record_from_solution = recordSolSECO,
            # record_from_solution = (x, p) -> norm(x[1:length(x)÷2] .- p, 8),
            )

sol = @time BifurcationKit.solve(prob, Newton(), NewtonPar(verbose = true, tol=1e-10))

plotsol(sol.u)

eig = EigArnoldiMethod(; sigma = 0.1, which = BK.LM(), tol = 1e-10)
opts = ContinuationPar(dsmin = 0.001, dsmax = 0.009, ds = -0.001, p_max = 3.65, p_min= 3.0, nev = 30, plot_every_step = 10, newton_options = NewtonPar( verbose = false, tol = 1e-11, eigsolver = eig, max_iterations = 10), max_steps = 100, n_inversion = 8, detect_bifurcation = 3)
br = @time continuation(prob, PALC(), opts;
        plot = true,
        verbosity = 2,
        normC = norminf,
    )

BK.plot(br)[1]

eigenvals(br, 1)
BifurcationKit.eigenvalsfrombif(br, 7)

bp = get_normal_form(br, 1; scaleζ = norminf, start_with_eigen = Val(false))
####################################################################################################
# automatic branch switching from bp point
T22 = [continuation(br, ind,
                ContinuationPar(BK.getcontparams(br); p_max = 3.7, p_min = 3.1, ds = -0.001, dsmax = 0.03, detect_bifurcation = 2, max_steps = 100);
                alg = PALC(
                            tangent = Bordered(), 
                            bls = BorderingBLS(;check_precision = true),
                            # dotθ = BifurcationKit.DotTheta((x,y) -> dot(x,y)/length(y)/2),
                            dotθ = BifurcationKit.DotTheta(DOT)

                ),
                start_with_eigen = Val(false),
                # callbackN  = BK.cbMaxNormAndΔp(0.1, 0.01),
                scaleζ = norminf,
                verbosity = 2,
                normC = norminf,
                plot = true) for ind in 1:10 if br.specialpoint[ind].type==:bp]

BK.plot(br, T22..., vars=(:param, :n2))[1]

T1 = continuation(br, 5,
                ContinuationPar(BK.getcontparams(br); p_max = 3.6, ds = -0.001, dsmax = 0.03, detect_bifurcation = 3, max_steps = 140, n_inversion = 6);
                alg = PALC(
                            # tangent = Bordered(),
                            bls = BorderingBLS(;check_precision = true),
                            # dotθ = BifurcationKit.DotTheta((x,y) -> dot(x,y)/length(y)/2),
                            dotθ = BifurcationKit.DotTheta(DOT)
                        ),
                start_with_eigen = Val(false),
                callbackN  = BifurcationKit.cbMaxNormAndΔp(0.1, 0.01),
                scaleζ = norminf,
                verbosity = 2,
                normC = norminf,
                # bothside = true,
                plot = true)

BK.plot(br, T22..., T1, vars=(:param, :n2), putspecialptlegend = false, branchlabel = fill("", 100))[1]
BK.plot(br, T1, putspecialptlegend = false, branchlabel = fill("", 100))[1]

opts_tss = ContinuationPar(BK.getcontparams(br), ds = -0.001, p_max = 3.8, dsmax = 0.025, dsmin=1e-4, detect_bifurcation = 3, max_steps = 500, plot_every_step = 50)
TSS = continuation(T1, 1,
                opts_tss;
                # alg = Natural(),
                # alg = MoorePenrose(tangent=PALC(tangent=Bordered()),
                alg = PALC(
                        tangent = Bordered(), 
                        bls = BorderingBLS(;check_precision = true),
                        # dotθ = BifurcationKit.DotTheta((x,y) -> dot(x,y)/length(y)/3),
                        dotθ = BifurcationKit.DotTheta(DOT)
                        ),
                verbosity = 2,
                normC = norminf,
                start_with_eigen = Val(false),
                callbackN  = BifurcationKit.cbMaxNormAndΔp(0.1, 0.01),
                plot = true)

BK.plot(TSS)[1]
BK.plot(br, T1, TSS, putspecialptlegend = false, branchlabel = fill("", 100))[1]
BK.plot(br, T22..., T1, TSS, vars=(:param, :n2), putspecialptlegend = false, branchlabel = fill("", 100))
####################################################################################################
opt_po = BK.NewtonPar(tol = 1e-8, max_iterations = 15)
opts_po_cont = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.001,
        p_max = 4.0, max_steps = 250, newton_options = opt_po,
        detect_bifurcation = 2, nev = 20, tol_stability = 1e-4, plot_every_step = 1)
prob_po = BK.Trapeze(; M = 40,
        jacobian = BK.FullSparseInplace(),
        # massmatrix = BK.getmassmatrix(prob, nothing, nothing)
        )
br_po = BK.continuation(TSS, 6, opts_po_cont, prob_po;
        alg = PALC(tangent = Bordered(), ),
        verbosity = 3,
        plot = true,
        start_with_eigen = Val(false),
        δp = 0.01,
        plot_solution = (ax, x, p; ax1, iter, state, k...) -> begin
            _sol = BK.get_periodic_orbit(BK.getprob(iter), x, p.p)
            heatmap!(ax, _sol.u'[:, 1:div(size(_sol.u,1),2)]; colormap = :viridis, k...)
        end,
        normC = BK.norminf)

BK.plot(br_po)[1]

####################################################################################################
# we could user a much simpler aBS from Hopf but tried a preconditioned solver
using FFTW

begin
    M = 40
    nf_hopf = BK.get_normal_form(br, 1; start_with_eigen = Val(false))
    pred = predictor(nf_hopf, .001; ampfactor = .91)
    list_of_time_steps = [pred.orbit(t) for t in LinRange(0, 2pi, M + 1)[1:M]]
    orbitguess_f = vcat(reduce(vcat, list_of_time_steps), pred.period)
    r_hopf = nf_hopf.params.j0
end

par_hopf = @set par_seco.j0 = r_hopf + 0.01

prob_po = BK.Trapeze(; M,
        prob_vf = BK.re_make(prob, params = par_hopf),
        N = length(prob.u0),
        xπ = zero(orbitguess_f[1:end-1]),
        ϕ = zero(orbitguess_f[1:end-1]),
        jacobian = BK.FullSparseInplace(),
        # massmatrix = BK.getmassmatrix(prob, nothing, nothing)
        )
BK.updatesection!(prob_po, orbitguess_f, par_hopf)

Prec = BK.POTrapCirculantPrec(prob_po, orbitguess_f, (par_hopf), ref=:average)
# ls_po = KrylovLSInplace(verbose = 0, rtol = 1e-3, atol = 1e-12, memory = 40, itmax = 50, Pl = Prec, ldiv = true, is_inplace = false, n = length(orbitguess_f), m = length(orbitguess_f))
ls_po = GMRESIterativeSolvers(verbose = false, reltol = 1e-3, restart = 40, maxiter = 50, Pl = Prec, log=true)

# ls_po = DefaultLS()
opt_po = NewtonPar(tol = 1e-7, linsolver = ls_po, verbose = true)
outpo_f = @time newton(prob_po, orbitguess_f, opt_po; normN = norminf)
Prec = BK.POTrapCirculantPrec(prob_po, orbitguess_f, (par_hopf), ref=:1)

opts_po_cont = BK.ContinuationPar(dsmin = 0.001, dsmax = 0.04, ds = 0.001,
        p_max = 4.0, max_steps = 250, newton_options = opt_po, p_min = 3.,
        detect_bifurcation = 0, nev = 20, tol_stability = 1e-4, plot_every_step = 10)

function record_po(x, p; iter, state, k...)
    (nrm = NormL2T(BK.getprob(iter), x, ustar(p.p.j0)),)
end

br_po = BK.continuation(br, 1, opts_po_cont, prob_po;
    alg = PALC(tangent = Secant(), bls = MatrixBLS()),
    verbosity = 3,
    # linear_algo = BorderingBLS(solver = ls_po, check_precision = false),
    linear_algo = MatrixFreeBLS(ls_po, false),
    plot = true,

    start_with_eigen = Val(false),
    δp = 0.01, ampfactor = 0.91,

    record_from_solution = record_po,
    plot_solution = (ax, x, p; ax1, iter, state, k...) -> begin
        _sol = BK.get_periodic_orbit(BK.getprob(iter), x, p.p)
        l1,l2 = size(_sol.u)
        # heatmap!(ax, _sol.t, 1:div(l1,2), _sol.u'[:, 1:div(l1,2)]; colormap = :viridis, k...)
        lines!(ax, _sol.t, _sol.u'[:, div(l1,4)])
    end,
    normC = BK.norminf)