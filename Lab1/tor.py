
import gmsh
import math
import os
import sys

gmsh.initialize()

gmsh.model.add("tor")

tor1 = gmsh.model.occ.addTorus(0, 0, 0, 3, 1, 7, 0.9*math.pi)
tor11 = gmsh.model.occ.addTorus(0, 0, 0, 3, 0.5, 1, 0.9*math.pi)
fluid1 = gmsh.model.occ.cut([(3, 7)], [(3, 1)])
gmsh.model.occ.synchronize()
tor1sem = gmsh.model.occ.symmetrize([(3, 7)], 0, -1, 0, 0)
gmsh.model.occ.synchronize()

tor2 = gmsh.model.occ.addTorus(0, 0, 0, 3, 1, 8, math.pi)
tor21 = gmsh.model.occ.addTorus(0, 0, 0, 3, 0.5, 2, math.pi)
fluid2 = gmsh.model.occ.cut([(3, 8)], [(3, 2)])
gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)

gmsh.write("tor.msh")
gmsh.write("tor.geo_unrolled")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()