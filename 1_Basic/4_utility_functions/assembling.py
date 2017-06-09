
from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = H1(mesh,order=3, dirichlet="top")

u, v = fes.TrialFunction(), fes.TestFunction()
bfi = SymbolicBFI(grad(u)*grad(v))
lfi = SymbolicLFI(v)

# use taskmanager to have parallel IterateElements
with TaskManager():
    u = MyAssemble(fes,bfi,lfi)

Draw(u)
