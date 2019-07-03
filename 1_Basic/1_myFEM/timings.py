from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import *
from myngspy import *
import time

SetNumThreads(4)
ngsglobals.msg_level=0

geometry = CSGeometry()
s = 0.3
outer_box = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1)).mat("outer")
inner_box = OrthoBrick(Pnt(-s,-s,-s),Pnt(s,s,s)).mat("inner")
geometry.Add(outer_box-inner_box)
geometry.Add(inner_box)

mesh = Mesh(geometry.GenerateMesh(maxh=0.2))

# for i in range(3):
#     mesh.Refine()

order = 10

cf = MyCoefficient()
Draw(cf, mesh, 'cf')
val = MyIntegrate(cf, mesh, order=order)
print('value', val)

# cflist = [cf if mat=='inner' else 0.0 for mat in mesh.GetMaterials()]
# cf = CoefficientFunction(cflist)
