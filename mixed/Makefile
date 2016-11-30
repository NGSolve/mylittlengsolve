objects = mixed_diffusion.o hybridmixed.o

%.o : %.cpp
	ngscxx -c $? -o $@

libmixed.so : $(objects)
	ngsld -shared $(objects) -lngcomp -lngsolve -o $@

clean:
	rm *.o libmixed.so

