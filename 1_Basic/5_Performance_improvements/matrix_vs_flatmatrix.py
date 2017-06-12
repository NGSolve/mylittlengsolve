from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *
from time import time

ngsglobals.msg_level = 1

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

fes = FESpace("myfespace", mesh, dirichlet="top|bottom|right|left")

tries = 1000

start = time()
with TaskManager():
    for i in range(tries):
        a = BilinearForm(fes)
        a += MyLaplace(CoefficientFunction(1))
        a.Assemble()
time_a = (time()-start)/tries

start = time()
with TaskManager():
    for i in range(tries):
        a2 = BilinearForm(fes)
        a2 += MyBetterLaplace(CoefficientFunction(1))
        a2.Assemble()
time_a2 = (time()-start)/tries

print("MyLaplace needed ", time_a*1000, " milliseconds for assembling")
print("MyBetterLaplace needed ", time_a2*1000, "milliseconds for assembling")
