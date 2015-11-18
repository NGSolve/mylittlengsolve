from ngsolve import *

from ctypes import CDLL
# on Windows replace '.so' with '.dll'
mylngs = CDLL("libmyngsolve.so")

m = Mesh("square.vol")

fes = FESpace("myfespace", m, dirichlet=[1,2,3,4], flags = { "secondorder" : True } )
print ("freedofs: ", fes.FreeDofs())

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += BFI("mylaplace", coef=1)
# a += SymbolicBFI(grad(u)*grad(v))

f = LinearForm(fes)
f += LFI("mysource", coef=1)
# f += SymbolicLFI(x*y*v)

a.Assemble()
f.Assemble()

u = GridFunction(fes)

print ("solve")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
