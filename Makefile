objects = all_in_one.o demo_instat.o demo_stokes.o myElement.o \
myHOElement.o myIntegrator.o demo_coupling.o demo_nonlinear.o  \
myFESpace.o myHOFESpace.o myAssembling.o


%.o : %.cpp
	gcc -O2 -fpic -DNETGEN_ELTRANS -I. -I$(NETGENDIR)/../include -c $? -o $@

libmyngsolve.so : $(objects)
	gcc -shared -fpic $(objects) -o $@

clean:
	rm *.o libmyngsolve.so



dist:
	cd ..; tar -czf MyLittleNGSolve-1.0.tar.gz my_little_ngsolve/Makefile my_little_ngsolve/*.cpp my_little_ngsolve/*.hpp my_little_ngsolve/*.pde my_little_ngsolve/*.vol my_little_ngsolve/*.in2d my_little_ngsolve/windows
