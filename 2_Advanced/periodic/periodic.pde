shared = libperiodic

mesh = periodic.vol.gz

define constant geometryorder = 3
define coefficient coef_source 0, 1

define fespace v -type=periodic_h1ho -order=3 -dirichlet=[1]

define gridfunction u -fespace=v

define bilinearform a -fespace=v -symmetric -fespace=v
laplace 1

define linearform f -fespace=v
source coef_source

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -solver=direct