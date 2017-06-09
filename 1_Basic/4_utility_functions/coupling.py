
from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes1 = H1(mesh, order=3)

u1,v1 = fes1.TrialFunction(), fes1.TestFunction()
a1 = BilinearForm(fes1)
a1 += SymbolicBFI(grad(u1)*grad(v1))
# penalty term
a1 += SymbolicBFI(1e7*u1*v1,definedon=mesh.Boundaries("top|right"))

f1 = LinearForm(fes1)
f1 += SymbolicLFI(v1)

u1 = GridFunction(fes1,"u1")

with TaskManager():
    a1.Assemble()
    f1.Assemble()
    u1.vec.data = a1.mat.Inverse() * f1.vec

fes2 = H1(mesh, order=6)
u2,v2 = fes2.TrialFunction(), fes2.TestFunction()
a2 = BilinearForm(fes2)
a2 += SymbolicBFI(grad(u2)*grad(v2))
# penalty term
a2 += SymbolicBFI(1e7*u2*v2,definedon=mesh.Boundaries("top|right"))

f2 = LinearForm(fes2)

# allocate the vector
f2.Assemble()

u2 = GridFunction(fes2,"u2")

with TaskManager():
    a2.Assemble()
    MyCoupling(u1,f2)
    u2.vec.data = a2.mat.Inverse() * f2.vec

Draw(u1)
Draw(u2)
