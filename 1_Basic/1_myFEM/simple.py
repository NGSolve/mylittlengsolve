from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

# for i in range(5):
#     mesh.Refine()

fes = MyFESpace(mesh, dirichlet="top|bottom|right|left") #, secondorder=True)
u,v = fes.TnT()
# print ("freedofs: ", fes.FreeDofs())


a = BilinearForm(fes)
# a += grad(u) * grad(v) * dx
a += MyLaplace(1)

f = LinearForm(fes)
# f += MyCoefficient()*v * dx
f += MySource(x*y)

u = GridFunction(fes)

# Exercise: Implement functionality for inhomogeneous boundary conditions
# This needs our finite element space to have a differential operator
# on the boundary.

with TaskManager(int(1e8)):
    a.Assemble()
    f.Assemble()
    print ("solve")
    u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
