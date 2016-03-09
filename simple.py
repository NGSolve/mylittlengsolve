from ngsolve import *

from ctypes import CDLL
# on Windows replace '.so' with '.dll'
mylngs = CDLL("libmyngsolve.so")

m = Mesh("square.vol")

# fes = FESpace("myfespace", m, dirichlet=[1,2,3,4], flags = { "secondorder" : True } )
fes = FESpace("myhofespace", m, dirichlet=[1,2,3,4], order = 5)
print ("freedofs: ", fes.FreeDofs())

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
# a += BFI("mylaplace", coef=1)
a += SymbolicBFI(grad(u)*grad(v))

f = LinearForm(fes)
# f += LFI("mysource", coef=x*y)
f += SymbolicLFI(x*y*v)

a.Assemble()
f.Assemble()

u = GridFunction(fes)

print ("solve")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
