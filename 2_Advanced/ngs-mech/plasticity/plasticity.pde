geometry = probe.in2d
mesh = probe.vol

define constant heapsize = 100000000

define constant geometryorder = 3

shared = ngs-mech


define coefficient one
1,


define coefficient coefE
1,

define coefficient coefnu
1,

define coefficient coefyield
1,

define coefficient coefg
0, 0.1, 0, 

define coefficient coefpen
1e6, 0, 0, 


define fespace v -type=plasticity -order=2

define gridfunction u -fespace=v

define bilinearform a -fespace=v -symmetric
plasticity coefE coefnu coefyield
robin coefpen -comp=1
robin coefpen -comp=2

define linearform f -fespace=v  -print
neumann coefg -comp=2



numproc plasticitysimple np1 -maxit=100



define bilinearform ux -nonassemble -fespace=v
mass one -comp=1

define bilinearform uy -nonassemble -fespace=v
mass one -comp=2

numproc drawflux dfux -bilinearform=ux -solution=u -label=ux
numproc drawflux dfuy -bilinearform=uy -solution=u -label=uy


define bilinearform p1 -nonassemble -fespace=v
mass one -comp=3

define bilinearform p2 -nonassemble -fespace=v
mass one -comp=4

define bilinearform p3 -nonassemble -fespace=v
mass one -comp=5


numproc drawflux dfp1 -bilinearform=p1 -solution=u -label=p1
numproc drawflux dfp2 -bilinearform=p2 -solution=u -label=p2
numproc drawflux dfp3 -bilinearform=p3 -solution=u -label=p3

numproc drawflux dfp4 -bilinearform=a -solution=u -label=sigma


define bilinearform bfdisp -nonassemble -fespace=v
plasticity_idu  

numproc drawflux dfp5 -bilinearform=bfdisp -solution=u -label=disp




