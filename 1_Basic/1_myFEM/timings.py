from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import *
from myngspy import *
import time

SetNumThreads(4)

geometry = CSGeometry()
s = 0.3
outer_box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).mat("outer")
inner_box = OrthoBrick(Pnt(-s,-s,-s),Pnt(s,s,s)).mat("inner")
geometry.Add(outer_box-inner_box)
geometry.Add(inner_box)

mesh = Mesh(geometry.GenerateMesh(maxh=0.2))

# for i in range(3):
#     mesh.Refine()

# Draw(mesh)
order = 10

def doTiming(mesh, cf):
    with TaskManager(pajetrace=10**8):
        t0 = time.time()
        val = MyIntegrate(cf, mesh, order=order)
        val2 = MyIntegrateSIMD(cf, mesh, order=order)
        val3 = MyIntegrateParallel(cf, mesh, order=order)
        val4 = Integrate(cf, mesh, order=order)
        print(val,val2, val3, val4)
        t1 = time.time()
    return t1-t0


cf = MyCoefficient()
cf = CoefficientFunction([cf if mat=='inner' else 0.0 for mat in mesh.GetMaterials()])
print('myCoefficient\t', doTiming(mesh,cf))
Draw(cf, mesh, 'MyCoefficient')
# print('x*y          \t', doTiming(mesh,x*y))
# print('x*y (comp)   \t', doTiming(mesh,(x*y).Compile()))
# print('x*y (comp true)\t', doTiming(mesh,(x*y).Compile(True,wait=True)))
