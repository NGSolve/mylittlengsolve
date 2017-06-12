from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = FESpace("myfespace", mesh, dirichlet="top|bottom|right|left")
# fes = FESpace("myfespace", mesh, dirichlet="top|bottom|right|left", flags={"secondorder":True})
print ("freedofs: ", fes.FreeDofs())

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += MyLaplace(CoefficientFunction(1))
# a += SymbolicBFI(grad(u)*grad(v))

f = LinearForm(fes)
f += MySource(x*y)
# f += SymbolicLFI(x*y*v)

a.Assemble()
f.Assemble()

u = GridFunction(fes)

print ("solve")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
