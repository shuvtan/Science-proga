import gmsh
import math
import os
import sys

gmsh.initialize()

gmsh.model.add("torus")
tor1 = gmsh.model.occ.addTorus(0, 0, 0, 3, 1, 1, 1.9*math.pi)
tor11 = gmsh.model.occ.addTorus(0, 0, 0, 3, 0.5, 11, 1.9*math.pi)
fluid = gmsh.model.occ.cut([(3, 1)], [(3, 11)])

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)

gmsh.write("torus.msh")
gmsh.write("torus.geo_unrolled")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()