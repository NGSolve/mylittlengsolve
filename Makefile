objects = all_in_one.o demo_instat.o demo_stokes.o myElement.o \
myHOElement.o myIntegrator.o demo_coupling.o demo_coupling_adv.o demo_nonlinear.o  \
myFESpace.o myHOFESpace.o myPreconditioner.o myAssembling.o linhypDG.o 

# linhypDG_par.o 


%.o : %.cpp
	ngscxx -I. -c $? -o $@


libmyngsolve.so : $(objects)
	ngscxx -shared $(objects) -lngsolve -o $@

clean:
	rm *.o libmyngsolve.so



dist:
	cd ..; tar -czf MyLittleNGSolve-5.1.tar.gz my_little_ngsolve/Makefile my_little_ngsolve/*.cpp my_little_ngsolve/*.hpp my_little_ngsolve/*.pde my_little_ngsolve/*.vol my_little_ngsolve/*.in2d my_little_ngsolve/windows/my_little_ngsolve* \
												  my_little_ngsolve/mixed/*.cpp my_little_ngsolve/mixed/*.hpp my_little_ngsolve/mixed*.pde \
												  my_little_ngsolve/mixed*.vol my_little_ngsolve/mixed*.in2d \
												  my_little_ngsolve/precond/*.cpp my_little_ngsolve/precond/*.hpp my_little_ngsolve/precond*.pde \
												  my_little_ngsolve/precond*.vol my_little_ngsolve/precond*.in2d \
												  my_little_ngsolve/tcl/*.cpp my_little_ngsolve/tcl/*.hpp my_little_ngsolve/tcl*.pde \
												  my_little_ngsolve/tcl*.vol my_little_ngsolve/tcl*.in2d \
												  my_little_ngsolve/windows/*.sln my_little_ngsolve/windows/*.vcproj
