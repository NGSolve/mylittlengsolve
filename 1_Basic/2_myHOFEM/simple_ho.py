from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = MyHOFESpace(mesh, dirichlet="top|bottom|right|left", order = 5)
u,v = fes.TnT()
print ("freedofs: ", fes.FreeDofs())

a = BilinearForm(fes)
a += grad(u)*grad(v)*dx

f = LinearForm(fes)
f += x*y*v*dx

a.Assemble()
f.Assemble()

u = GridFunction(fes)

print ("solve")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
