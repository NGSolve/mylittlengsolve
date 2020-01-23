geometry = pendel_constraint.geo
mesh = pendel_constraint.vol


shared = mbs

define coefficient e
1e0, 
2.1e3, 2.1e3, 2.1e3,

define coefficient nu
0.3, 0.3, 0.3,

define coefficient rho
7.8e-6, 7.8e-6, 7.8e-6,

define coefficient penalty
1e8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

#define coefficient coef_source
#(-0.05), 0, 0,

define coefficient coef_load
(-1e-2), (0), (0)



define fespace v -dim=3 -order=2
define fespace vscal -dim=1 -order=2
define fespace vp -order=1 -dim=6

define gridfunction u -fespace=v
define gridfunction u0 -fespace=v
define gridfunction usmall -fespace=v
define gridfunction v -fespace=v
define gridfunction veca -fespace=v
define gridfunction stress -fespace=vp

define bilinearform a -fespace=v -symmetric
elasticity e nu


define bilinearform m -fespace=v -symmetric
mass rho -comp=1
mass rho -comp=2
mass rho -comp=3



define bilinearform hmat -fespace=v -symmetric

define linearform f -fespace=v
source coef_load -comp=3



define preconditioner c -type=multigrid -bilinearform=a


numproc mbs np1 -bilinearforma=a -bilinearformm=m  -gridfunctionu=u -gridfunctionv=v -gridfunctionusmall=usmall -gridfunctionstress=stress -dt=2e-3 -te=20 -connections=[1,0] 

