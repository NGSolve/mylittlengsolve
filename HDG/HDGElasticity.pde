# geometry = beam.in2d
# mesh = beam.vol.gz
 geometry = plate.geo
 mesh = plate.vol.gz

# define constant testout = "test.out"
# define constant numthreads = 1
define constant heapsize = 20000000

shared = libHDG

define constant E = 1000
define constant nu = 0.3
define coefficient cf  (0,0,1)

define fespace v -type=HDGElasticity -order=3 -dirichlet=[1]  

define gridfunction u -fespace=v -addcoef

define linearform f -fespace=v
sourceedge cf -comp=1

define bilinearform a -fespace=v -symmetric -noprintelmat -xxelmatev
HDG_elasticity E nu

numproc bvp np1 -gridfunction=u -bilinearform=a -linearform=f -solver=direct

