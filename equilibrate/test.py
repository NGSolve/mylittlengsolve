from ngsolve import *
from netgen.geom2d import unit_square
from equilibrate import Equilibrate

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
ngsglobals.testout = "test.out"

fes = H1(mesh, order=3, dirichlet=[1,2,3,4])
u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()

source = CoefficientFunction(1)
f = LinearForm(fes)
f += SymbolicLFI(source*v)
f.Assemble()

u = GridFunction(fes, name="solution")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

primal_flux = grad(u)
eqflux = Equilibrate(primal_flux, source, fes, order=3)

Draw(u)
Draw(primal_flux, mesh, "primal_flux")
Draw(eqflux)
Draw(div(eqflux), mesh, "div_eqflux")

