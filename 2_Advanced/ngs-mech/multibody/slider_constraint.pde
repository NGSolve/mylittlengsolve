geometry = slider_constraint.geo
mesh = slider_constraint.vol


shared = mbs

# define constant geometryorder = 2

define coefficient e
2.1e5, 2.1e5, 2.1e5,

define coefficient nu
0.3, 0.3, 0.3,

define coefficient rho
7.8e-6, 7.8e-6, 7.8e-6,

define coefficient penalty
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

#define coefficient coef_source
#(-0.05), 0, 0,

define coefficient coef_load
(-0.2), (0), (0)

define coefficient one
1, 1, 1

define fespace v -dim=3 -order=1 -hb
define fespace vscal -dim=1 -order=1 -hb
define fespace vp -dim=6 -order=1

define gridfunction u -fespace=v
define gridfunction u0 -fespace=v
define gridfunction usmall -fespace=v
define gridfunction v -fespace=v
define gridfunction veca -fespace=v
define gridfunction stress -fespace=vp

define bilinearform a -fespace=v
elasticity e nu


define bilinearform m -fespace=v
mass rho -comp=1
mass rho -comp=2
mass rho -comp=3

define bilinearform hmat -fespace=v

define linearform f -fespace=v
source coef_load -comp=2

define bilinearform aid -fespace=v -nonassemble
mass one

numproc drawflux npdraw1 -bilinearform=a -solution=usmall -applyd -label=stress2

# numproc setvisual npvis -mminval=-1000000 -mmaxval=1000000 -autoscale=0 -vecfunction=u -deformation=1

numproc mbs np1 -bilinearforma=a -bilinearformm=m  -gridfunctionu=u -gridfunctionv=v -gridfunctionusmall=usmall -gridfunctionstress=stress -dt=0.005 -te=2 -TI=R2 -connections=[1,0,2,0,3,4,5,6] -connectionsx=[7,0,8,0,9,10,12,0] -connectionsy=[7,0,8,0,12,0]


