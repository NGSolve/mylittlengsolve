ngs-mech.so: plasticity_simple.cpp
	gcc -O2 -shared -fpic  -DNETGEN_ELTRANS -DLAPACK -I. -I/opt/netgen/include  $? -o $@
