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

function plotgridap!(ax, x::Vector; k...) 
    uh = zero(U)
    uh.free_values .= x
    plotgridap!(ax, uh; k...)
end

plotgridap!(ax, x; k...) = plot!(ax, Ω, x)

####################################################################################################
function buildDisk(h::Float64, R, filename = "MyDisk.msh")
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
buildDisk(N::Int, R, filename = "MyDisk.msh") = buildDisk(√(2*R*pi/N), R, filename)

buildDisk(800, 1)
####################################################################################################
using GridapGmsh
model = GmshDiscreteModel("MyDisk.msh")
const Ω = Triangulation(model)
num_cells(Ω)
GLMakie.plot(Ω)
GLMakie.wireframe(Ω, color=:black, linewidth=2)
GLMakie.scatter(Ω, marker=:star8, markersize=20, color=:blue)
#############################################
# function spaces
begin
order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, conformity=:H1,)#dirichlet_tags="boundary")
U = TrialFESpace(V)

# Ω = Triangulation(model)
degree = 2*order
dΩ = Measure(Ω, degree)

NL(u) = exp(u)
res(u, p, v) = ∫( -∇(v)⋅∇(u) -  v ⋅ (u - p.λ ⋅ (NL ∘ u)) * 10 )*dΩ
jac(u, p, du, v) = ∫( -∇(v)⋅∇(du) - v ⋅ du ⋅ (1 - p.λ *( NL ∘ u)) * 10 )*dΩ
d2res(u, p, du1, du2, v) = ∫( v ⋅ du1 ⋅ du2 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ
d3res(u, p, du1, du2, du3, v) = ∫( v ⋅ du1 ⋅ du2 ⋅ du3 ⋅ (NL ∘ u) * 10 * p.λ )*dΩ

uh = zero(U)
par_bratu = (λ = 0.01,)

# weight for normbratu
const w = cumsum(ones(length(uh.free_values))) / length(uh.free_values)
n = length(uh.free_values)
w .= (1 .+ LinRange(-1,1,n)) |> vec
w .-= minimum(w)
normbratu(x) = norm(x .* w) / sqrt(length(x))

prob = GridapBifProblem(res, uh, par_bratu, V, U, (@optic _.λ);
                jac,
                d2res = d2res,
                d3res = d3res,
                plot_solution = (ax,x,p; k...) -> plotgridap!(ax, x;  k...),
                record_from_solution = (x, p; k...) -> normbratu(x))
end
# factorize leads pivots issues, better use LU factorization here
optn = NewtonPar(eigsolver = EigArpack())#EigKrylovKit(dim = 100))
sol = BifurcationKit.solve(prob, Newton(), NewtonPar(optn; verbose = true))

opts = ContinuationPar(p_max = 40., p_min = 0.01, ds = 0.01, max_steps = 1000, detect_bifurcation = 3, newton_options = optn, nev = 20, tol_stability = 1e-6, n_inversion = 6)
br = continuation(prob, PALC(tangent = Bordered()), opts;
    plot = true,
    verbosity = 1,
    )

BifurcationKit.plot(br)

nf = get_normal_form(br, 4; verbose = true, scaleζ = norminf, autodiff = false)
####################################################################################################
br1 = continuation(br, 4,
        ContinuationPar(opts; ds = 0.001, dsmax = 0.05, max_steps = 140, detect_bifurcation = 3);
        verbosity = 1, plot = true, nev = 10,
        # usedeflation = true,
        scaleζ = norminf,
        autodiff = false,
        callback_newton = BifurcationKit.cbMaxNorm(100),
        )

BifurcationKit.plot(br, br1...)
####################################################################################################
diagram = @time bifurcationdiagram(prob, PALC(),
    # important argument: this is the maximal
    # recursion level
    3,
    ContinuationPar(opts; ds = 0.001, dsmax = 0.05, max_steps = 140, detect_bifurcation = 3);
    verbosity = 0, plot = true,
    # callback_newton = cb,
    # usedeflation = true,
    # finalise_solution = finSol,
    autodiff = false,
    verbosediagram = true,
    normC = norminf)

fig, ax = BifurcationKit.plot(diagram)
BifurcationKit.plot!(ax, br1)
BifurcationKit.plot!(ax, br1[2])
fig