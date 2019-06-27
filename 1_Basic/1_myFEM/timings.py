from netgen.csg import unit_cube
from ngsolve import *
from myngspy import *
import time

mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))

# for i in range(3):
#     mesh.Refine()

# Draw(mesh)

def doTiming(mesh, cf):
    with TaskManager(pajetrace=10**8):
        t0 = time.time()
        val = MyIntegrate(cf, mesh, order=2)
        val2 = Integrate(cf, mesh, order=2)
        print(val,val2)
        t1 = time.time()
    return t1-t0


print('myCoefficient\t', doTiming(mesh,MyCoefficient()))
print('x*y          \t', doTiming(mesh,x*y))
print('x*y (comp)   \t', doTiming(mesh,(x*y).Compile()))
print('x*y (comp true)\t', doTiming(mesh,(x*y).Compile(True,wait=True)))
