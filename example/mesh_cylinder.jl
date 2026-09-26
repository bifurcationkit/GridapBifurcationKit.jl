# example/mesh_cylinder.jl
#
# Mesh for the 2-D compressible flow past a cylinder example of
#   Douglas & Jolivet, Comput. Phys. Commun. 326 (2026) 110221 (Section 3.3)
# mirroring the ff-bifbox builder
#   examples/douglas_jolivet_2026/mesh_cylinder.md
#
# Geometry (all lengths in units of the cylinder diameter):
#   rectangle [-100, 150] x [0, 100]  minus  the half disk of radius 0.5 at the
#   origin. The full (mirrored) mesh uses [-100, 150] x [-100, 100] minus the
#   full disk.
#
# Boundary physical groups:
#   "inlet"  (x = xmin), "out" (x = xmax), "wall" (cylinder),
#   "lat"    (y = +rmax, and y = -rmax for the full mesh),
#   "sym"    (y = 0, half mesh only).
#
# The target element sizes reproduce the FreeFEM border discretizations:
#   cylinder   ~ 1/n0   (= 0.02 for n0 = 50)
#   symmetry   ~ 1/n1   (= 0.5  for n1 = 2)
#   far field  ~ 1/n2   (= 4.0  for n2 = 0.25)
#
# Run:  julia --project=example example/mesh_cylinder.jl

using GridapGmsh, gmsh_jll

const gmsh = GridapGmsh.gmsh

# Domaine full (cylinder_full.msh, utilisé par CylinderGridap.jl)
#  y
#  ↑                                                       centre = origine (0,0)
# 100 ┌───────────────────────────────────────────────────────────────┐ ← tag "lat"
#     │                                                               │
#     │            champ lointain   (taille ≈ size_far = 1/n2)        │
#     │                                                               │
#  20 ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· · ← bsds ≠ 0
#     │       couche absorbante  bsds = 1e-4[max(x-50,-20-x,0)²      │
#     │                                              + max(|y|-20,0)²] │
#     │                 ┌───────────────────────────┐                │
#   0 │─────────────────│   bsds = 0                │────────────────│ ← y=0 (sillage)
#     │        ●        │   x∈[-20,50], y∈[-20,20]  │                │
#     │     cylindre    └───────────────────────────┘                │
#     │     rayon rcyl = 0.5   ↑ raffiné le long de y=0 (MathEval)   │
#     │     tag "wall", taille ≈ size_cyl = 1/n0                     │
# -20 ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·
#     │                                                               │
# -100└───────────────────────────────────────────────────────────────┘ ← tag "lat"
#    xmin = -100                     x = 0                    xmax = 150
#    tag "inlet"               (100 amont)              (150 aval)  tag "out"
# Demi-domaine (cylinder.msh / cylinder_test.msh)
#  y
#  ↑
# 100 ┌──────────────────────────────────────────────┐ ← tag "lat"
#     │              champ lointain                   │
#  20 ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·· ·
#     │                                              │
#   0 │─────────────●────────────────────────────────│ ← tag "sym"
#     │          cylindre (r = 0.5, "wall")          │
#    (pas de y < 0)
#    xmin = -100                              xmax = 150
#    "inlet"                                   "out"

# Classify a boundary curve from its bounding box.
function _classify_curve(bb, xmin, xmax, rmax)
    x0, y0, _, x1, y1, _ = bb
    if x0 > -0.6 && x1 < 0.6 && y0 > -0.6 && y1 < 0.6
        return :wall
    elseif abs(x0 - xmin) < 1e-6 && abs(x1 - xmin) < 1e-6
        return :inlet
    elseif abs(x0 - xmax) < 1e-6 && abs(x1 - xmax) < 1e-6
        return :out
    elseif y0 > rmax - 1e-6 || y1 < -rmax + 1e-6
        return :lat
    elseif y0 > -1e-6 && y1 < 1e-6
        return :sym
    else
        return :other
    end
end

"""
    write_cylinder_msh(; full=false, n0=50, n1=2, n2=0.25, xmin=-100, xmax=150,
                       rmax=100, rcyl=0.5, dir=joinpath(@__DIR__, "meshes"),
                       name=full ? "cylinder_full" : "cylinder")

Write a triangular `.msh` (Gmsh 4.1) for the half (`full=false`) or mirrored
full (`full=true`) cylinder domain, with the physical groups listed in the file
header. Returns the file path.
"""
function write_cylinder_msh(; full = false, n0 = 50, n1 = 2, n2 = 0.25,
                            xmin = -100.0, xmax = 150.0, rmax = 100.0, rcyl = 0.5,
                            dir = joinpath(@__DIR__, "meshes"),
                            name = full ? "cylinder_full" : "cylinder")
    isdir(dir) || mkpath(dir)
    file = joinpath(dir, name * ".msh")
    size_cyl = 1.0 / n0
    size_axis = 1.0 / n1
    size_far = 1.0 / n2

    gmsh.clear()
    gmsh.model.add(name)
    occ = gmsh.model.occ

    ylo = full ? -rmax : 0.0
    hgt = full ? 2 * rmax : rmax
    rect = occ.addRectangle(xmin, ylo, 0.0, xmax - xmin, hgt)
    disk = occ.addDisk(0.0, 0.0, 0.0, rcyl, rcyl)
    out, _ = occ.cut([(2, rect)], [(2, disk)])
    occ.synchronize()
    surf = out[1][2]

    # --- boundary classification and physical groups -----------------------
    groups = Dict{Symbol, Vector{Int}}()
    for (dim, tag) in gmsh.model.getBoundary([(2, surf)], false, false, false)
        bb = gmsh.model.getBoundingBox(dim, tag)
        key = _classify_curve(bb, xmin, xmax, rmax)
        key === :other && continue
        push!(get!(groups, key, Int[]), tag)
    end
    for key in (:inlet, :out, :wall, :lat, :sym)
        haskey(groups, key) || continue
        gmsh.model.addPhysicalGroup(1, groups[key], -1, string(key))
    end
    gmsh.model.addPhysicalGroup(2, [surf], -1, "domain")

    # --- element size fields ----------------------------------------------
    fm = gmsh.model.mesh.field
    d1 = fm.add("Distance")
    fm.setNumbers(d1, "CurvesList", Float64.(groups[:wall]))
    fm.setNumber(d1, "Sampling", 200)
    t1 = fm.add("Threshold")
    fm.setNumber(t1, "InField", d1)
    fm.setNumber(t1, "SizeMin", size_cyl)
    fm.setNumber(t1, "SizeMax", size_far)
    fm.setNumber(t1, "DistMin", 0.5)
    fm.setNumber(t1, "DistMax", 25.0)
    fields = [t1]

    if haskey(groups, :sym)
        d2 = fm.add("Distance")
        fm.setNumbers(d2, "CurvesList", Float64.(groups[:sym]))
        fm.setNumber(d2, "Sampling", 200)
        t2 = fm.add("Threshold")
        fm.setNumber(t2, "InField", d2)
        fm.setNumber(t2, "SizeMin", size_axis)
        fm.setNumber(t2, "SizeMax", size_far)
        fm.setNumber(t2, "DistMin", 1.0)
        fm.setNumber(t2, "DistMax", 10.0)
        push!(fields, t2)
    else
        # full mesh: refine near the wake centerline y = 0
        m = fm.add("MathEval")
        fm.setString(m, "F",
            "min($(size_axis) + ($(size_far) - $(size_axis))*abs(y)/10.0, $(size_far))")
        push!(fields, m)
    end

    fmin = fm.add("Min")
    fm.setNumbers(fmin, "FieldsList", Float64.(fields))
    fm.setAsBackgroundMesh(fmin)

    # --- mesh and write ----------------------------------------------------
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.model.mesh.generate(2)
    gmsh.write(file)
    return file
end

"""
    generate_cylinder_meshes(; dir=joinpath(@__DIR__, "meshes"), kwargs...)

Generate the half (`cylinder.msh`) and full (`cylinder_full.msh`) meshes.
"""
function generate_cylinder_meshes(; dir = joinpath(@__DIR__, "meshes"), kwargs...)
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        return [write_cylinder_msh(; full = false, dir, kwargs...),
                write_cylinder_msh(; full = true, dir, kwargs...)]
    finally
        gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    for f in generate_cylinder_meshes()
        println("wrote ", f)
    end
end
