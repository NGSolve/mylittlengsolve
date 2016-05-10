NGSCXX = ${NETGENDIR}/ngscxx
NGSLD = ${NETGENDIR}/ngsld


all: demo_std demo_bla  demo_fem  demo_comp  demo_solve demo_parallel
# demo_comp1d

demo_std:  demo_std.cpp
	$(NGSCXX) -c demo_std.cpp 
	$(NGSLD) -o $@ demo_std.o -lngstd

demo_bla:  demo_bla.cpp
	$(NGSCXX) -c demo_bla.cpp 
	$(NGSLD) -o $@ demo_bla.o -lngstd -lngbla -lblas

demo_fem:  demo_fem.cpp
	$(NGSCXX) -c demo_fem.cpp 
	$(NGSLD) -o $@ demo_fem.o -lngfem -lngstd -lngbla

demo_comp:  demo_comp.cpp
	$(NGSCXX) -c demo_comp.cpp 
	$(NGSLD) -o $@ demo_comp.o -lngcomp -lngsolve -lngla -lngfem  -lngstd -lnglib -linterface

demo_comp1d:  demo_comp1d.cpp netgen_dummy.cpp
	$(NGSCXX) -c demo_comp1d.cpp netgen_dummy.cpp
	$(NGSLD) -o $@ demo_comp1d.o netgen_dummy.o -o demo_comp1d -lngcomp -lngsolve -lngla -lngfem  -lngstd -lnglib -linterface




demo_solve:  demo_solve.cpp
	$(NGSCXX) -c demo_solve.cpp 
	$(NGSLD) -o $@ demo_solve.o -lngcomp -lngsolve -lngla -lngfem  -lngstd -lnglib -linterface


demo_parallel:  demo_parallel.cpp
	$(NGSCXX) -c demo_parallel.cpp 
	$(NGSLD) -o $@ demo_parallel.o -lngstd


install:

clean:
	rm demo_std demo_bla demo_fem demo_comp demo_solve demo_parallel
