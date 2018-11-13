from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = MyFESpace(mesh, secondorder=True)

u = GridFunction(fes,"shapes")

Draw(u)

# we can use the additionally exported function here
for i in range(fes.ndof):
    print("Draw basis function ", i)
    u.vec[:] = 0
    u.vec[i] = 1
    Redraw()
    input("press key to draw next shape function")
