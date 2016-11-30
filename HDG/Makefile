objects = HDGLaplace.o HDGElasticity.o

%.o : %.cpp
	ngscxx -c $? -o $@

libHDG.so : $(objects)
	ngsld -shared $(objects) -lngcomp -lngsolve -o $@

clean:
	rm *.o libHDG.so

