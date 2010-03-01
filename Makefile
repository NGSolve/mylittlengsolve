objects = all_in_one.o demo_instat.o demo_stokes.o myElement.o \
myHOElement.o myIntegrator.o demo_coupling.o demo_nonlinear.o  \
myFESpace.o myHOFESpace.o


%.o : %.cpp
	gcc -O2 -fpic -DNETGEN_ELTRANS -I. -I$(NETGENDIR)/../include -c $? -o $@

libmyngsolve.so : $(objects)
	gcc -shared -fpic $(objects) -o $@

clean:
	rm *.o libmyngsolve.so
