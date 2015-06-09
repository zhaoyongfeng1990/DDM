ddm : functions.o main.o ISFSwimmer.o ISFRunAndTumbleAndDiffusion.o ISFRunAndTumble2D.o ISFRunAndTumble3D.o ISFRTDWithP.o ISFSwimmerSimp.o
	g++-5 *.o -o ddm -O3 -lm -lfftw3 -lgslcblas -lgsl -fopenmp
functions.o : functions.cpp functions.h
	g++-5 -O3 -c functions.cpp -std=gnu++11
main.o : main.cpp functions.h
	g++-5 -O3 -c main.cpp -std=gnu++11
ISFSwimmer.o : ISFSwimmer.cpp functions.h
	g++-5 -O3 -c ISFSwimmer.cpp -std=gnu++11
ISFRunAndTumbleAndDiffusion.o : ISFRunAndTumbleAndDiffusion.cpp functions.h
	g++-5 -O3 -c ISFRunAndTumbleAndDiffusion.cpp -std=gnu++11
ISFRunAndTumble2D.o : ISFRunAndTumble2D.cpp functions.h
	g++-5 -O3 -c ISFRunAndTumble2D.cpp -std=gnu++11
ISFRunAndTumble3D.o : ISFRunAndTumble3D.cpp functions.h
	g++-5 -O3 -c ISFRunAndTumble3D.cpp -std=gnu++11
ISFRTDWithP.o : ISFRTDWithP.cpp functions.h
	g++-5 -O3 -c ISFRTDWithP.cpp -std=gnu++11
ISFSwimmerSimp.o : ISFSwimmerSimp.cpp functions.h
	g++-5 -O3 -c ISFSwimmerSimp.cpp -std=gnu++11
debugmain.o : debugmain.cpp functions.h
	g++-5 -O3 -c debugmain.cpp -std=gnu++11
clean :
	rm ./*.o
recover :
	make ddm
	./ddm recover
debug : debugmain.o functions.o ISFSwimmer.o ISFRunAndTumbleAndDiffusion.o ISFRunAndTumble2D.o ISFRunAndTumble3D.o ISFRTDWithP.o
	g++-5 debugmain.o functions.o ISFSwimmer.o ISFRunAndTumbleAndDiffusion.o ISFRunAndTumble2D.o ISFRunAndTumble3D.o ISFRTDWithP.o -o debug -O3 -lm -lfftw3 -lgslcblas -lgsl -fopenmp
	./debug
