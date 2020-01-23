geometry = balken.geo
mesh = balken.vol

define constant localheap = 10000000

shared = libngsmech

define constant geometryorder = 1

define coefficient e
1e11,

define coefficient e2
2e11,

define coefficient nu

0.2,

define coefficient coef_force_y
0, 3e9, 0,


define fespace v -h1ho -order=3 -dim=3 -dirichlet=[1]
define gridfunction u -fespace=v -nested
define gridfunction ulin -fespace=v -nested

define fespace vstress -h1ho -order=3 -dim=6
define gridfunction stress -fespace=vstress

define bilinearform alin -fespace=v -symmetric
elasticity e2 nu

define bilinearform a -fespace=v -symmetric
nonlinmech e2 nu

define linearform f -fespace=v
neumann coef_force_y -comp=2
# neumann coef_force_x -comp=1

define preconditioner c -type=direct -bilinearform=alin

numproc bvp np1lin -bilinearform=alin -linearform=f -gridfunction=ulin -preconditioner=c 

numproc nonlinelast np1 -bilinearform=a -linearform=f -gridfunction=u -maxsteps=50


numproc drawflux np2 -bilinearform=a -solution=u -label=elstress
numproc calcflux np3 -bilinearform=a -solution=u -flux=stress
