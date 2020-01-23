.cpp.o:
	ngscxx -O2 -c -I.  $? -o $@


largestrain.so: nonlinelast.o  nonlinmech_integrator.o
	ngsld -shared -fpic nonlinelast.o nonlinmech_integrator.o -o libngsmech.so -L/opt/netgen/lib -lngsolve

