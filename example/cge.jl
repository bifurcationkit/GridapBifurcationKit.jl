cd(@__DIR__)
using Pkg
pkg"activate ."

using Revise
# using Plots
using GridapGmsh
using GridapGmsh: gmsh
using GridapMakie, GLMakie
Makie.inline!(true)
using Gridap
using Gridap.FESpaces
using GridapBifurcationKit
using BifurcationKit
const BK = BifurcationKit
####################################################################################################
function buildDisk(h::Float64, R, filename = "MyDisk-CGE.msh")
        gmsh.initialize()
        model = gmsh.model
    
        # Parameters
        # R = 1    # Radius
        # h = 0.021 # Mesh size
        center = model.geo.addPoint(0, 0, 0, h, -1)
    
        # Create 3 Points on the circle
        points = []
        Nd = 3
        for j in 1:Nd
            push!(points, model.geo.addPoint(R*cos(2*pi*j/Nd), R*sin(2*pi*j/Nd), 0, h))
        end
    
        # Create 3 circle arc
        lines = []
        for j in 1:Nd
            push!(lines, model.geo.addCircleArc(points[j], center, points[j+1 < Nd+1 ? j + 1 : 1]))
        end
    
        # Curveloop and Surface
        curveloop = model.geo.addCurveLoop([1,2,3])
        disk = model.geo.addPlaneSurface([curveloop])
    
        # This command is mandatory and synchronize CAD with GMSH Model. The less you launch it, the better it is for performance purpose
        gmsh.model.geo.synchronize()
    
        # Physical groups
        # gmsh.model.addPhysicalGroup(dim, list of tags, physical tag)
        t1 = gmsh.model.addPhysicalGroup(1, lines, 1)
        gmsh.model.addPhysicalGroup(2, [disk], 10)
    
        gmsh.model.setPhysicalName(1, t1, "boundary")
    
        # Mesh (2D)
        model.mesh.generate(2)
        # Write on disk
        gmsh.write(filename)
    
        # run the Gui
        # gmsh.fltk.run();
        # Finalize GMSH
        gmsh.finalize()
end

# function  where we specify a rough number of triangles
buildDisk(N::Int, R, filename = "MyDisk-CGE.msh") = buildDisk(√(2*R*pi/N), R, filename)

buildDisk(250, pi, "MyDisk-CGE.msh")
####################################################################################################
using GridapGmsh
model = GmshDiscreteModel("MyDisk-CGE.msh")
const Ω = Triangulation(model)
num_cells(Ω)
GLMakie.plot(Ω)
GLMakie.wireframe(Ω, color=:black, linewidth=2)
GLMakie.scatter(Ω, marker=:star8, markersize=6, color=:blue)

begin
    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe, conformity=:H1,)  # homogeneous natural (Neumann) BC, like the pde2path cge demo
    U = TrialFESpace(V)

    # 2 fields: (u1, u2) = (Re A, Im A)
    Y = MultiFieldFESpace([V, V])
    const X = MultiFieldFESpace([U, U])

    # Ω = Triangulation(model)
    degree = 2*order
    dΩ = Measure(Ω, degree)
end
#############################################
function plotgridap!(ax, x::Vector; k...) 
    uh = zero(X)
    uh.free_values .= x
    u1, _ = uh
    plotgridap!(ax, u1; k...)
end

function plotgridap!(ls, x; k...)
    hm=plot!(ls, Ω, x)
    gc = ls.layoutobservables.gridcontent[]
    gp = gc === nothing ? ls.figure[1, 2] :
         gc.parent[gc.span.rows, maximum(gc.span.cols) + 1]
    Colorbar(gp, hm)
end

plotgridap(x; k...) = (f = Figure(); plotgridap!(Axis(f[1,1]), x; k...); f)
####################################################################################################
# Complex Ginzburg-Landau equation (CGE), cubic-quintic:
#
#   ∂_t A = r A + i δ ν A + Δ A - (c3 + i μ)|A|² A - c5 |A|⁴ A
#
# with A = u1 + i u2. The comoving advection `s*(-y ∂x + x ∂y)` enters the
# residual:
#   f1 = r u1 - ν u2       - ua(c3 u1 - μ u2) - c5 ua² u1
#   f2 = r u2 + δ ν u1     - ua(c3 u2 + μ u1) - c5 ua² u2
# where ua = |A|² and ∂θ = -y ∂x + x ∂y.
#
# The residual `res` below equals -G = M f - K u - s Krot u (weak form), so that
# BifurcationKit's convention M ∂_t u = res coincides with the pde2path
# convention M ∂_t u = -G. `jac` is the corresponding derivative of `res`.

# advection field of the comoving frame: (-y, x) = ∂θ
adv(x) = VectorValue(-x[2], x[1])

function res((u1,u2), p, (v1,v2))
    ua(U1, U2) = U1^2 + U2^2
    f1(U1, U2) = p.r*U1 -      p.ν *U2 - ua(U1,U2)*(p.c3*U1 - p.μ*U2) - p.c5*ua(U1,U2)^2*U1
    f2(U1, U2) = p.r*U2 + (p.δ*p.ν)*U1 - ua(U1,U2)*(p.c3*U2 + p.μ*U1) - p.c5*ua(U1,U2)^2*U2
    ∫(((           f1∘(u1,u2))*v1 
                + (f2∘(u1,u2))*v2 
                - ∇(u1)⋅∇(v1) 
                - ∇(u2)⋅∇(v2)
                - p.s*((adv⋅∇(u1))*v1 + (adv⋅∇(u2))*v2) ))*dΩ
end

function jac((u1,u2), p, (du1,du2), (v1,v2))
    ua(X1, X2) = X1^2 + X2^2
    d1f1(X1, X2) = p.r - 2X1*(p.c3*X1 - p.μ*X2) - p.c3*ua(X1,X2) - p.c5*(ua(X1,X2)^2 + 4X1^2*ua(X1,X2))
    d2f1(X1, X2) = -p.ν - 2p.c3*X1*X2 + 2p.μ*X2^2 + p.μ*ua(X1,X2) - 4p.c5*ua(X1,X2)*X1*X2
    d1f2(X1, X2) = p.δ*p.ν - 2p.c3*X1*X2 - 2p.μ*X1^2 - p.μ*ua(X1,X2) - 4p.c5*ua(X1,X2)*X1*X2
    d2f2(X1, X2) = p.r - 2X2*(p.c3*X2 + p.μ*X1) - p.c3*ua(X1,X2) - p.c5*(ua(X1,X2)^2 + 4X2^2*ua(X1,X2))
    ∫( ((d1f1∘(u1,u2))*du1*v1 + (d2f1∘(u1,u2))*du2*v1 
                - ∇(du1)⋅∇(v1)
                + (d1f2∘(u1,u2))*du1*v2 + (d2f2∘(u1,u2))*du2*v2 
                - ∇(du2)⋅∇(v2)
                - p.s*((adv⋅∇(du1))*v1 + (adv⋅∇(du2))*v2)) )*dΩ
end

generateur1((u1,u2), p, (du1,du2), (v1,v2)) = ∫( adv⋅∇(du1)*v1 + adv⋅∇(du2)*v2) * dΩ
generateur2((u1,u2), p, (du1,du2), (v1,v2)) = ∫(  -du2*v1 + du1*v2  )*dΩ

function N2(x)
    u1 = zero(X)
    u1.free_values .= x
    n = sum(∫(u1[1]*u1[1]+u1[2]*u1[2])dΩ)
    s = sum(∫(1)dΩ)
    sqrt(n/s)
end

function N∞(x; R = pi*0.9, npts = 300)
    uh = zero(X)
    uh.free_values .= x
    u1, u2 = uh
    pts = [Gridap.Point(xx, yy) for xx in LinRange(-R, R, npts),
                           yy in LinRange(-R, R, npts) if xx^2 + yy^2 ≤ R^2]
    maximum(sqrt.(evaluate(u1, pts).^2 .+ evaluate(u2, pts).^2))
end

amplitude(x) = maximum(x) - minimum(x)

# mass matrix (identity on (u1, u2), nonsingular)
m0((u1,u2), (v1,v2)) =  ∫(u1*v1 + u2*v2)*dΩ

# parameters (inspired by the pde2path csh demo): r continued from -0.1
par_cge = (r = -0.1, ν = 1.0, μ = 5., c3 = -1., c5 = 1.0, δ = 1.0, s = 0.0)

# Lie generators
Gθ     = Gridap.FESpaces.assemble_matrix((du, v) -> generateur1(zero(X), par_cge, du, v), X, Y)
Gphase = Gridap.FESpaces.assemble_matrix((du, v) -> generateur2(zero(X), par_cge, du, v), X, Y)

uh = zero(X); uh.free_values .= 0.01#rand(length(uh.free_values))
prob = GridapBifProblem(res, uh, par_cge, Y, X, dΩ, (@optic _.r);
                        jac = jac, 
                        # d2res = d2res, 
                        # d3res = d3res, 
                        mass = m0,
                        plot_solution = (ax, x, p; ax1, k...) -> plotgridap!(ax, x; k...),
                        record_from_solution = (x, p; k...) -> (N = N2(x),)
                        )

eig_gridap = EigArnoldiMethod(sigma = 0.4, which  = BK.LM())
eig_gridap = EigArpack(sigma = 0.4, which  = :LM)
optn = NewtonPar(eigsolver = eig_gridap)
sol = @time BK.solve(prob, BK.Newton(), NewtonPar(optn; verbose = true))
plotgridap(sol.u )

opts = ContinuationPar(p_max = 0.5, p_min = -0.5, ds = 0.005, dsmax = 0.05, max_steps = 1000,
                       detect_bifurcation = 3, nev = 30, n_inversion = 8, newton_options = optn,
                       plot_every_step = 2)
br = @time BK.continuation(prob, BK.PALC(), opts; verbosity = 0, plot = true)
BK.plot(br)[1]
####################################################################################################
# plotgridap(br.eig[6].eigenvecs[:,4]|>real)
# nf = BK.get_normal_form(br, 1; start_with_eigen = Val(false), verbose = true)
# plotgridap(imag(nf.ζ))
####################################################################################################
# SW: A=B
# RW  A > B = 0
function guessFromHopfO2(branch, ind_hopf, eigsolver, M, A, B;  k = 1, phase = 0)
    specialpoint = branch.specialpoint[ind_hopf]
    ind_ev = specialpoint.ind_ev
    @error "" ind_ev branch.eig[specialpoint.idx].eigenvals
    # parameter value at the Hopf point
    p_hopf = specialpoint.param
    # frequency at the Hopf point
    ωH = abs(imag(branch.eig[specialpoint.idx].eigenvals[ind_ev]))/k
    # eigenvectors for the eigenvalues iω
    ζ0 = geteigenvector(eigsolver, br.eig[specialpoint.idx][2], ind_ev)
    ζ0 ./=  norm(ζ0)
    ζ1 = geteigenvector(eigsolver, br.eig[specialpoint.idx][2], ind_ev - 2)
    ζ1 ./=  norm(ζ1)
    @error "" dot(Gθ*ζ0, ζ1) dot(Gθ*ζ0, ζ0) dot(Gphase*ζ0, ζ1) dot(Gphase*ζ0, ζ0)
    orbitguess = [real.(specialpoint.x .+
                 A .* ζ0 .* cis(t - phase) .+
                 B .* ζ1 .* cis(t - phase)) for t in LinRange(0,2pi,M+1)[1:end]]
    return (; p = p_hopf, period = 2pi/ωH, guess = orbitguess, x0 = specialpoint.x, ζ0 = ζ0, ζ1 = ζ1, ω = ωH)
end

M = 10 # number of time slices (plotting purposes)
pred = guessFromHopfO2(br, 2, optn.eigsolver, M, 40*(1+im), 40*(1+im); k = 2); #TW
r_hopf = pred.p; Th = pred.period; orbitguess_tw = pred.guess; ωH = pred.ω
u_ref = copy(orbitguess_tw[1][1:length(prob.u0)])

begin
    f = Figure()
    n = 1
    for i=1:3,j=1:3
        ga = f[i, j] = GridLayout()
        plotgridap!(Axis(ga[1,1]), real.(orbitguess_tw[n][1:length(prob.u0)]) )
        n += 1
    end
    f
end

probTW = TWModel(re_make(prob, params = (par_cge..., r = r_hopf - 0.01)), 
            # (Gθ, Gphase), 
            Gθ,
            # Gphase,
            u_ref; 
            DAE = [false],
            jacobian = BK.FullLU())

plotgridap(probTW.u₀)

begin
    wave = newton(probTW, 
            # vcat(u_ref, 0),
            vcat(u_ref, ωH),
            # vcat(u_ref, 0, ωH),
            NewtonPar(verbose = true, max_iterations = 25, tol = 1e-9),
            normN = norminf,
        )
    println("s = ", wave.u[end-1], ", ω = " , wave.u[end], ", ||u||₂ = ", N2(wave.u[1:end-probTW.nc]))
    plotgridap(wave.u[1:end-probTW.nc] )
end

opt_cont_br = ContinuationPar(p_min = -5.5, p_max = 1.5, newton_options = optn, ds= 0.001, dsmax = 0.04, detect_bifurcation = 3, nev = 30, n_inversion = 8, tol_stability = 1e-5, plot_every_step = 2, max_steps = 100)
br_TW = @time continuation((probTW), deepcopy(wave.u), PALC(tangent = Bordered()), opt_cont_br;
    eigsolver = BK.EigenWave(br.contparams.newton_options.eigsolver, false),
    record_from_solution = (x, p; k...) -> (max = maximum(x[1:end-probTW.nc]), nrm = norm(x[1:end-probTW.nc]), s = x[end], amp = amplitude(x[1:end-probTW.nc])),
    plot_solution = (ax, x, p; ax1, k...) -> begin
                    plotgridap!(ax, x[1:end-probTW.nc]; k...)
                    # BK.plot!(ax1, br)
    end,
    finalise_solution = (z, tau, step, contResult; k...) -> begin
            N2(z.u[1:end-probTW.nc]) > 0.001
            true
    end, 
    verbosity = 1,
    normC = norminf,
    plot = true,
    # bothside = true,
    )

begin
    f, ax = BK.plot(br_TW, vars = (:param, :))
    BK.plot!(ax, br)
    f
end

# - #  1,     hopf at r ≈ +0.60878789 ∈ (+0.60878120, +0.60878789), |δp|=7e-06, [converged], δ = (-2, -2), step =  14
# - #  2,     hopf at r ≈ +0.77953149 ∈ (+0.77953064, +0.77953149), |δp|=8e-07, [converged], δ = (-2, -2), step =  18
# - #  3,     hopf at r ≈ +1.24308802 ∈ (+1.24297898, +1.24308802), |δp|=1e-04, [converged], δ = ( 2,  2), step =  27
######################################################################
using SparseArrays

hp = BK.get_normal_form(br_TW, 3; start_with_eigen = Val(false), detailed = Val(true), verbose = true)

probMTW = Trapeze(M = 20 ;
    N = length(br_TW.prob.u0),
    # massmatrix = spdiagm(0 => vcat(ones(length(prob.u0)), 0.)),
    # linear solver for the periodic orbit problem
    # OPTIONAL, one could use the default
    jacobian = BK.FullLU())

opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds= 0.005, p_max = 2.05, max_steps = 20, newton_options = optn, nev = 30, tol_stability = 1e-2, detect_bifurcation = 2, plot_every_step = 1)
opts_po_cont = @set opts_po_cont.newton_options.max_iterations = 15
opts_po_cont = @set opts_po_cont.newton_options.tol = 1e-5
opts_po_cont = @set opts_po_cont.newton_options.verbose = true
 
br_MTW1 =  continuation(
    br_TW, 2,
    opts_po_cont,
    probMTW ;
    eigsolver = BK.FloquetQaD(DefaultEig(), false),
    verbosity = 3,

    δp = 0.002, #ampfactor = 0.9,
    start_with_eigen = Val(false),
    # use_normal_form = false, 
    
    # tangent predictor
    alg = Natural(),
    # regular parameters for the continuation
    # a few parameters saved during run
    record_from_solution = (u, p; k...) -> begin
        outt = BK.get_periodic_orbit(p.prob, u, p)
        m = maximum(outt.u[eachindex(uh.free_values),:])
        nrm = norm(outt.u[eachindex(uh.free_values),:])
        return (;max = m, nrm , s = u[end-1], period = u[end])
    end,
    # plotting of a section
    plot_solution = (ax, x, p; ax1, k...) -> begin
        outt = BK.get_periodic_orbit(p.prob, x, p)
        lines!(ax, outt.t, outt.u[end, :])
        BK.plot!(ax1, br_TW, vars = (:param, :max))
    end,
    # print the Floquet exponent
    finalise_solution = (z, tau, step, contResult; k...) -> begin
        true
    end,
    plot = true,
    callback_newton = BK.cbMaxNorm(2),
    normC = norminf
    )

BK.plot(br_TW, br_MTW1)

br_MTW1.period

lines(br_MTW1.param, abs.(2pi ./ br_MTW1.s) ./ br_MTW1.period)

BK.plot(br_MTW1)[1]

begin
    ind = 20
    sol = BK.get_periodic_orbit(br_MTW2, ind)
    f = Figure(size=(700,700))
    inds = round.(Int, LinRange(1, size(sol.u, 2), 9))   # 9 tranches temporelles
    for (k, (i, j)) in enumerate(Iterators.product(1:3, 1:3))
        ga = f[i, j] = GridLayout()
        plotgridap!(Axis(ga[1,1]), real.(sol.u[1:length(prob.u0),inds[k]]))
    end
    f
end


br_MTW2 =  continuation(
    br_TW, 3,
    opts_po_cont,
    probMTW ;
     eigsolver = BK.FloquetQaD(DefaultEig(), false),
    verbosity = 3,

    δp = 0.002, #ampfactor = 0.9,
    start_with_eigen = Val(false),
    # use_normal_form = false, 
    
    # tangent predictor
    alg = PALC(),
    # regular parameters for the continuation
    # a few parameters saved during run
    record_from_solution = (u, p; k...) -> begin
        outt = BK.get_periodic_orbit(p.prob, u, p)
        m = maximum(outt.u[eachindex(uh.free_values),:])
        nrm = norm(outt.u[eachindex(uh.free_values),:])
        return (;max = m, nrm , s = u[end-1], period = u[end])
    end,
    # plotting of a section
    plot_solution = (ax, x, p; ax1, k...) -> begin
        outt = BK.get_periodic_orbit(p.prob, x, p)
        lines!(ax, outt.t, outt.u[end, :])
        BK.plot!(ax1, br_TW, vars = (:param, :max))
    end,
    # print the Floquet exponent
    finalise_solution = (z, tau, step, contResult; k...) -> begin
        true
    end,
    plot = true,
    callback_newton = BK.cbMaxNorm(2),
    normC = norminf
    )

BK.plot(br_TW, br_MTW2)
br_MTW2.period
lines(br_MTW2.param, abs.(2pi ./ br_MTW2.s) ./ br_MTW2.period)
BK.plot(br_MTW2)[1]
######################################################################
pred = guessFromHopfO2(br, 2, optn.eigsolver, M, 40*(1+0im), 0*(1+im); k = 1); #SW
r_hopf = pred.p; Th = pred.period; orbitguess_sw = pred.guess; ωH = pred.ω
u_ref = copy(orbitguess_sw[1][1:length(prob.u0)])

begin
    f = Figure()
    n = 1
    for i=1:3,j=1:3
        ga = f[i, j] = GridLayout()
        plotgridap!(Axis(ga[1,1]), real.(orbitguess_sw[n][1:length(prob.u0)]) )
        n += 1
    end
    f
end

probSW = TWModel(re_make(prob, params = (par_cge..., r = r_hopf - 0.0)), 
            # (Gθ, Gphase), 
            # Gθ,
            Gphase,
            deepcopy(u_ref); 
            DAE = [0],
            update_section_every_step = 1,
            jacobian = BK.FullLU())

plotgridap(probSW.u₀)

(probSW.∂u₀) |> norminf
        
# BK.residual_tw!(probSW, vcat(0*u_ref, 0.), vcat(u_ref, -ωH), probSW.prob_vf.params)|>lines

begin
    optn = NewtonPar(verbose = true, max_iterations = 40, tol = 1e-10)
    wave_sw = newton(probSW, 
        vcat(u_ref, 1*0.03),
        # vcat(u_ref, -ωH, 0ωH),
        # vcat(u_ref, 0, ωH),
        optn,
        normN = norminf,
    )
println("s = ", wave_sw.u[end-1], ", ω = " , wave_sw.u[end], ", ||u||₂ = ", N2(wave_sw.u[1:end-probSW.nc]))

f=Figure()
plotgridap!(Axis(f[1,1]), u_ref )
plotgridap!(Axis(f[1,3]), wave_sw.u[1:end-probSW.nc] )
f
end

amplitude(x) = maximum(x) - minimum(x)
opt_cont_br = ContinuationPar(p_min = -5.5, p_max = 1.5, newton_options = NewtonPar(optn, max_iterations = 5), ds= -0.001, dsmax = 0.05, detect_bifurcation = 3, nev = 10, n_inversion = 6, tol_stability = 1e-3, plot_every_step = 2, max_steps = 130)
br_SW = @time continuation((probSW), deepcopy(wave_sw.u), PALC(tangent = Bordered()), opt_cont_br;
    eigsolver = BK.EigenWave(br.contparams.newton_options.eigsolver, false),
    record_from_solution = (x, p; k...) -> (max = maximum(x[1:end-probSW.nc]), nrm = norm(x[1:end-probSW.nc]), s = x[end], amp = amplitude(x[1:end-probSW.nc])),
    plot_solution = (ax, x, p; ax1, k...) -> begin
                    plotgridap!(ax, x[1:end-probSW.nc]; k...)
                    # BK.plot!(ax1, br)
    end,
    finalise_solution = (z, tau, step, contResult; k...) -> begin
            N2(z.u[1:end-probSW.nc]) > 0.001
            true
    end, 
    verbosity = 2,
    normC = norminf,
    plot = true,
    # bothside = true,
    )

begin
    f, ax = BK.plot(br_SW, br_TW, vars = (:param, :max))
    BK.plot!(ax, br)
    f
end

BK.plot(br, br_SW, br_TW, br_MTW1, br_MTW2)