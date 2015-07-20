CXX = g++-5
CXXFLAGS = -fopenmp -O3 -std=c++11 -DHAVE_INLINE
LFLAGS = -lm -lgsl -lgslcblas -O3 -lfftw3 -lfftw3l -fopenmp -std=c++11

Objects = aveQ.o averSqrModTau.o clean.o ddm.o fitting.o loadFiles.o main.o print.o NILT.o NILT_RTDP.o ISFRTD.o ISFSWP.o ISFSW.o ISFRTDP.o
debugObjects = ISFRTD.o ISFSWP.o ISFSW.o ISFRTDP.o aveQ.o averSqrModTau.o clean.o ddm.o fitting.o loadFiles.o debugmain.o print.o NILT.o NILT_RTDP.o

ddm : $(Objects)
	$(CXX) -o ddm $(LFLAGS) $(Objects)

debug : $(debugObjects)
	$(CXX) -o debug $(LFLAGS) $(debugObjects)
	./debug

$(Objects) : ddm.h parameters.h ISFtype.h
$(debugObjects) : ddm.h parameters.h

debugmain.o NILT.o NILT_RTDP.o ISFRTD.o ISFRTDP_NoLT_sigma.o : NILT.h

.PHONY : clean
clean :
	-rm *.o