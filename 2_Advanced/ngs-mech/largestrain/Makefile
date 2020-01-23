.cpp.o:
	gcc -O2 -c -fpic -DUSE_TIMEOFDAY -I. -DLAPACK -I/opt/netgen/include  $? -o $@


largestrain.so: nonlinelast.o  nonlinmech_integrator.o
	gcc -shared -fpic nonlinelast.o nonlinmech_integrator.o -o libngsmech.so -L/opt/netgen/lib -lngsolve

