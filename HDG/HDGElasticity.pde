# geometry = beam.in2d
# mesh = beam.vol.gz
# geometry = plate.geo
# mesh = plate.vol.gz
geometry = pipe.geo
mesh = pipe.vol.gz

define constant heapsize = 20000000
define constant geometryorder = 2

shared = libHDG

define constant E = 1000
define constant nu = 0.2
define coefficient cf  (0,0,1)

define fespace v -type=HDGElasticity -order=2 -dirichlet=[1]  

define gridfunction u -fespace=v -addcoef

define coefficient cg 0, 1, 0, 
define linearform f -fespace=v
sourceedge cf -comp=1
# neumannhdiv cg -comp=2

define bilinearform a -fespace=v -symmetric -noprintelmat -xxelmatev
HDG_elasticity E nu

numproc bvp np1 -gridfunction=u -bilinearform=a -linearform=f -solver=direct

