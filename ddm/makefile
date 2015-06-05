ddm : main.cpp functions.cpp functions.h
	g++-5 *.o -o ddm -O3 -lm -lfftw3 -lgslcblas -lgsl -fopenmp
functions.o : functions.cpp functions.h
	g++-5 -O3 -c functions.cpp -std=gnu+11
main.o : main.cpp functions.h
	g++-5 -O3 -c main.cpp -std=gnu+11
ISFSwimmer.o : ISFSwimmer.cpp functions.h
	g++-5 -O3 -c ISFSwimmer.cpp -std=gnu+11
ISFRunAndTumbleAndDiffusion.o : ISFRunAndTumbleAndDiffusion.cpp functions.h
	g++-5 -O3 -c ISFRunAndTumbleAndDiffusion.cpp -std=gnu+11
ISFRunAndTumble2D.o : ISFRunAndTumble2D.cpp functions.h
	g++-5 -O3 -c ISFRunAndTumble2D.cpp -std=gnu+11
ISFRunAndTumble3D.o : ISFRunAndTumble3D.cpp functions.h
	g++-5 -O3 -c ISFRunAndTumble3D.cpp -std=gnu+11
clean :
	rm ./ddm
recover :
	make ddm
	./ddm recover