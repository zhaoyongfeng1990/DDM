CXX = g++-5
CXXFLAGS = -fopenmp -O3 -std=c++11 -I/usr/local/include -DHAVE_INLINE
LFLAGS = -L/usr/local/lib -lm -lgsl -lgslcblas -O3 -lfftw3 -lfftw3l -fopenmp -std=c++11

Objects = ddm.o fitting.o loadFiles.o main.o print.o NILT.o NILT_RTDP.o ISFRTDP.o ISFRTDPfix.o ISFRTDPTT.o ISFRTDPTTfix.o
debugObjects = ISFRTDP.o ISFRTDPfix.o ISFRTDPTT.o ISFRTDPTTfix.o ddm.o fitting.o loadFiles.o debugmain.o print.o NILT.o NILT_RTDP.o

ddm : $(Objects)
	$(CXX) -o ddm $(LFLAGS) $(Objects)

debug : $(debugObjects)
	$(CXX) -o debug $(LFLAGS) $(debugObjects)
	./debug

$(Objects) : ddm.h parameters.h ISFtype.h
$(debugObjects) : ddm.h parameters.h

debugmain.o NILT.o NILT_RTDP.o ISFRTDP.o ISFRTDPfix.o ISFRTDPTT.o ISFRTDPTTfix.o : NILT.h

.PHONY : clean
clean :
	-rm *.o
