from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = FESpace("myhofespace", mesh, order=5)

u = GridFunction(fes,"shapes")

Draw(u)

def printshape(i):
    print("Draw basis function ", i)
    u.vec[:] = 0
    u.vec[i] = 1
    Redraw()

print("ndof = ", fes.ndof)
for i in range(len(u.vec)):
    printshape(int(input("enter number of shapefunction to print:")))
