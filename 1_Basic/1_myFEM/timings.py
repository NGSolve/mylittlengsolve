from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import *
from myngspy import *
import time

SetNumThreads(1)
ngsglobals.msg_level=0

geometry = CSGeometry()
s = 0.3
outer_box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).mat("outer")
inner_box = OrthoBrick(Pnt(-s,-s,-s),Pnt(s,s,s)).mat("inner")
geometry.Add(outer_box-inner_box)
geometry.Add(inner_box)

mesh = Mesh(geometry.GenerateMesh(maxh=0.2))

for i in range(1):
    mesh.Refine()

order = 5

cf = MyCoefficient()
Draw(cf, mesh, 'cf')
val = MyIntegrate(cf, mesh, order=order)
val2 = Integrate(cf, mesh, order=order)
print('values', val, val2)

# cflist = [cf if mat=='inner' else 0.0 for mat in mesh.GetMaterials()]
# cf = CoefficientFunction(cflist)

for t in Timers():
    if 'Integrate' in t['name']:
        print(t['name'],t['counts'],t['time'])
