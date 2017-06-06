from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = MyFESpace(mesh, dirichlet="top|bottom|right|left", order = 5)
print ("freedofs: ", fes.FreeDofs())

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
# a += BFI("mylaplace", coef=1)
a += SymbolicBFI(grad(u)*grad(v))

f = LinearForm(fes)
# f += LFI("mysource", coef=x*y)
f += SymbolicLFI(x*y*v)

a.Assemble()
f.Assemble()

u = GridFunction(fes)

print ("solve")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
