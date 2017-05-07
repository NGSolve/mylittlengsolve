objects = myngspy.o

%.o : %.cpp
	ngscxx -I. -c $? -o $@


libmyngspy.so : $(objects)
	ngsld -v -shared $(objects) -lngfem -lngcomp -lngsolve -o $@

clean:
	rm *.o libmyngspy.so

