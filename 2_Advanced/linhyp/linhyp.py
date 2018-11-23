from ngsolve import *
from liblinhyp import *
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
fes = L2(mesh, order=5, all_dofs_together=True)

gfu = GridFunction(fes)
u0 = exp(-20**2 * ( (x-0.7)**2 + (y-0.5)**2 ))
gfu.Set (u0)
Draw (gfu, autoscale=False, min=0, max=1)

wind = CoefficientFunction( (y-0.5, 0.5-x) )

if True:
    # our own convection
    conv = Convection(fes, wind)
else:
    # convection from NGS-Py
    u,v = fes.TnT()
    bn = wind*specialcf.normal(mesh.dim)
    conv = BilinearForm(fes)
    conv += SymbolicBFI (u * wind*grad(v))
    conv += SymbolicBFI ( (-bn*IfPos(bn, u, u.Other()) * (v-v.Other())), VOL, skeleton=True)
    conv += SymbolicBFI ( (-bn*IfPos(bn, u, 0) * v), BND, skeleton=True)


t = 0
tau = 5e-3
tend = 10

w = gfu.vec.CreateVector()
hu = gfu.vec.CreateVector()

with TaskManager():
    while t < tend:
        # improved Euler's method
        conv.Apply (gfu.vec, w)
        fes.SolveM (rho=1, vec=w)
        hu.data = gfu.vec + 0.5*tau*w

        conv.Apply (hu, w)
        fes.SolveM (rho=1, vec=w)
        gfu.vec.data += tau * w
        t += tau
        Redraw()
    

