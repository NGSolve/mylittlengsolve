# Test the new exported coefficientfunction

from netgen.geom2d import unit_square
from ngsolve import *
from mylittlengs import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
cf = MyCoefficient("coefficientvalue.txt")
Draw(cf,mesh,"myCoefficient")
