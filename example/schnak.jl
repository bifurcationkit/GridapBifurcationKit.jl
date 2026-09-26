using Revise
using GLMakie, Gridap, SparseArrays, BifurcationKit, Arpack
using LinearAlgebra
using Gridap.FESpaces
using GridapBifurcationKit
Makie.inline!(true)
const BK = BifurcationKit

# plotsol(u, Nx, Ny) = heatmap(reshape(u[1:Nx*Ny], Nx, Ny))
# plotsol(u::SingleFieldFEFunction, Nx, Ny) = heatmap(reshape(map(x->x[3], u.cell_vals), Nx, Ny))
# plotsol!(u, Nx, Ny; k...) = heatmap!(reshape(u[1:Nx*Ny], Nx, Ny);k...)
# plotgridap!(x; k...) = (l=length(x)÷2;n=Int(sqrt(l));heatmap!(reshape(x[l+1:end], n, n); color=:viridis, k...))
# plotgridap(x; k...) =( plot();plotgridap!(x; k...))
#############################################
# discretisation
begin
Nx = 20
kc = √( √2 - 1)
lx = 2pi / kc
ly = lx / sqrt(3)
Ny = round(ly/lx*Nx) |> Int 
domain = (-lx, lx, -ly, ly)
# domain = (0, 2lx, 0, 2ly)
cells = (Nx, Ny)
model = CartesianDiscreteModel(domain, cells)

# using GridapGmsh
# model = GmshDiscreteModel("/Users/rveltz/Downloads/crisscross2.msh")

# function spaces
order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model,reffe,conformity=:H1,) #dirichlet_tags="boundary")
U = TrialFESpace(V)

Y = MultiFieldFESpace([V, V])
X = MultiFieldFESpace([U, U])

Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω, degree)
end

function res((u,v),p,(V1, V2))
    NL1(u,v) =  -u + u^2*v + p.σ * (u-1/v)^2
    NL2(u,v) = p.λ - u^2*v - p.σ * (u-1/v)^2

      ∫(  -∇(u)⋅∇(V1) + NL1∘(u,v) ⋅ V1 +
    (-p.d) * ∇(v)⋅∇(V2) + NL2∘(u,v) ⋅ V2 )*dΩ
end

function jac((u,v),p,(du, dv),(V1, V2))
    d1NL1(u,v) = -1 + 2*u*v + 2*p.σ*(u - 1/v)
    d2NL1(u,v) = u^2 + 2*p.σ*(u - 1/v)/v^2
    d1NL2(u,v) = -2*u*v - 2*p.σ*(u - 1/v)
    d2NL2(u,v) = -u^2 - 2*p.σ*(u - 1/v)/v^2

    ∫(  -∇(du)⋅∇(V1) - p.d * ∇(dv)⋅∇(V2) +
            d1NL1∘(u,v)⋅V1⋅du + d2NL1∘(u,v)⋅V1⋅dv +
            d1NL2∘(u,v)⋅V2⋅du + d2NL2∘(u,v)⋅V2⋅dv )*dΩ
end



par_sh = (λ = 3.35, σ = 0.0, d = 60.0)
uh = zero(X)
uh.free_values .= 1
####################################################################################################
function plotsol!(ax, sol; nd = 100, k...)
    solr, soli = sol

    lx1 = minimum(x->x[1], model.grid_topology.vertex_coordinates)
    lx2 = maximum(x->x[1], model.grid_topology.vertex_coordinates)

    ly1 = minimum(x->x[2], model.grid_topology.vertex_coordinates)
    ly2 = maximum(x->x[2], model.grid_topology.vertex_coordinates)

    X = LinRange(lx1,lx2,nd)
    Y = LinRange(ly1,ly2,nd)
    ur = [evaluate(solr, Gridap.Point(x,y)) for x in X, y in Y]
    ui = [evaluate(soli, Gridap.Point(x,y)) for x in X, y in Y]
    Makie.heatmap!(ax, X, Y, ur.^2 .+ui.^2; k...)
end
plotsol2!(ax, sol::Vector; k...) = (_uh = zero(X); _uh.free_values .= sol; plotsol!(ax, _uh; k...))
plotsol2(sol; k...) = (plot();plotsol2!(sol; k...))
plotsol(sol; k...) = (plot();plotsol!(sol; k...))
####################################################################################################
using Statistics

recordSolSCH(x, p; k...) = (
            # u0 = reshape(x[1:length(x)÷2],Nx,Ny)[Nx÷2,Ny÷2],
            u0 = x[length(x)÷4],
            nrm =  norm(x[1:length(x)÷2], Inf),
            max = maximum(x[1:length(x)÷2]),
            min = minimum(x[1:length(x)÷2]),
            zero = maximum(x[1:length(x)÷2]) - p,
            x = norm(x[1:length(x)÷2] .* w) - p * norm(w),
            norminf = norm(x[1:length(x)÷2], Inf))

prob = GridapBifProblem(res, uh, par_sh, Y, X, dΩ, (@optic _.λ);
            jac,
            plot_solution = (ax,x,p;ax1, iter, state, kp...) -> plotsol2!(ax, x; kp...),
            record_from_solution = recordSolSCH,
            )

const w = rand(length(BK.getu0(prob))÷2)

eig = EigArpack(0.1, :LM, tol = 1e-12)
# eig = EigArnoldiMethod(sigma = 0.1, which  = LM())
# eig = EigKrylovKit()
optn = NewtonPar(eigsolver = eig)

sol = @time BK.solve(prob, Newton(), NewtonPar(optn; verbose = true, tol=1e-11))

opts = ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = -0.01, p_max = 3.5, p_min= 2., detect_bifurcation = 3, nev = 30, newton_options = NewtonPar(optn; verbose = false, tol = 1e-11), max_steps = 100, tol_stability = 1e-7, n_inversion = 4)
br = @time continuation(prob, PALC(tangent = Bordered()), opts;
        plot = true,
        # verbosity = 2,
    )

BK.plot(br)[1]

eigenvals(br, 13) |> display
####################################################################################################
bp = get_normal_form(br, 1; verbose = true, nev = 20, scaleζ = norminf, start_with_eigen = Val(false))

br1 = @time continuation(br, 1, ContinuationPar(br.contparams; ds = 1e-2, max_steps = 100, detect_bifurcation = 3, dsmax = 0.02, dsmin = 1e-3, plot_every_step = 5, p_min = 2.2, n_inversion = 6);
        # verbosity = 3, 
        plot = true,
        scaleζ = norminf,
        start_with_eigen = Val(false),
        )

BK.plot(br, br1...; dash_unstable_style = true, vars = (:param, :x))[1]

begin
f,ax = BK.plot(br; vars = (:param, :max),label="")
cols = Makie.wong_colors()
for (i, b) in enumerate(br1)
    c = cols[mod1(i + 1, length(cols))]
    BK.plot!(ax, b; vars = (:param, :max), branchcolor = c, label="", )
    BK.plot!(ax, b; vars = (:param, :min), branchcolor = c, label="", )
end

xlims!(ax,(2.25,3.4))
# ylims!(ax,(1,4.25))
# legend_bottom_right!(f)
f
end
####################################################################################################
inds_br_hopf = [ii for ii in eachindex(br1) if any(sp -> sp.type == :hopf, br1[ii].specialpoint)]
brh = br1[inds_br_hopf[1]]
ind_h = findfirst(sp -> sp.type == :hopf, brh.specialpoint)

opt_po = BK.NewtonPar(tol = 1e-8, verbose = true, max_iterations = 15)
opts_po_cont = BK.ContinuationPar(dsmin = 0.001,
        dsmax = 0.04, ds = 0.01,
        p_max = 3.0,
        max_steps = 100,
        newton_options = opt_po,
        plot_every_step = 1,
        detect_bifurcation = 0,
        nev = 20,
        tol_stability = 1e-6,
        )

hp = BK.get_normal_form(brh, ind_h; start_with_eigen = Val(false))
# number of time slices for the periodic orbit
M = 41
probFD = BK.Trapeze(;M,
    # specific method for solving linear system
    # of Periodic orbits with trapeze method
    # You could use the default one :FullLU (slower here)
    jacobian = BK.FullSparseInplace(),
    N = length(sol.u),
    # massmatrix = BK.getmassmatrix(prob, nothing, nothing),
    )

pred = predictor(hp, .1; ampfactor = 0.1)
list_of_time_steps = [pred.orbit(t) for t in LinRange(0, 2pi, M + 1)[1:M]]
orbitguess_f = vcat(reduce(vcat, list_of_time_steps), pred.period)

begin
    f = Figure()
    ax = Axis(f[1,1])
    solt = BK.get_periodic_orbit(poTrap, orbitguess_f, nothing)
    lines!(ax, solt.t, solt.u[size(solt.u,1)÷4, :])
    f
end

poTrap = Trapeze(re_make(prob, params = (@set par_sh.λ = hp.p + 0.01)), 
                    BK.residual(prob, list_of_time_steps[1], (@set par_sh.λ = hp.p + 0.01)) |> normalize, 
                    zeros(length(prob.u0)), 
                    M, 
                    length(prob.u0);
                    # jacobian = BK.FullMatrixFree(), 
                    jacobian = BK.FullSparseInplace(),
                    )

using FFTW
Prec = BK.POTrapCirculantPrec(poTrap, orbitguess_f, (@set par_sh.λ = hp.p - 0.01), ref=:1);
# Prec = @time ilu(Jpo, τ = 0.005);

ls_po = GMRESIterativeSolvers(verbose = false, reltol = 1e-4, restart = 40, maxiter = 50, Pl = Prec, log=true)
# ls_po = KrylovLSInplace(;n=size(Jpo,1), m=size(Jpo,1), is_inplace=false, Pl = Prec, verbose = 0, ldiv = true, rtol = 1e-3);
# @reset opts_po_cont.newton_options.linsolver = ls_po

br_po = continuation(brh, ind_h,
            opts_po_cont,
            probFD
            ;
            start_with_eigen = Val(false),
            δp = 0.01,
            ampfactor = 0.02,
            verbosity = 2,
            plot = true,
)

####################################################################################################
diagram = @time bifurcationdiagram(prob, br,
    # recursion level:    2,
    (args...) -> ContinuationPar(opts; ds = 0.001, dsmax = 0.05, max_steps = 140, detect_bifurcation = 3);
    verbosity = 0, plot = true,
    usedeflation = true,
    normC = norminf)

f, ax = BK.plot(diagram)
xlims!(ax, (3,3.25));f
