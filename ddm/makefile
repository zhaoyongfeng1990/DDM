CXX = g++-5
CXXFLAGS = -fopenmp -O3 -std=c++11 -DHAVE_INLINE
Objects = ISFRTDWithP.o ISFRunAndTumble2D.o ISFRunAndTumble3D.o ISFRunAndTumbleAndDiffusion.o ISFRunAndTumbleAndDiffusionNoLT.o ISFSwimmer.o ISFSwimmerSimp.o LaplaceTrans.o aveQ.o averSqrModTau.o clean.o ddm.o fitting.o loadFiles.o main.o print.o shiftImage.o NILT.o NILT_RTDP.o ISFRTDP_NoLT.o ISFRTDP_NoLT_sigma.o
debugObjects = debugmain.o ISFRTDWithP.o ISFRunAndTumble2D.o ISFRunAndTumble3D.o ISFRunAndTumbleAndDiffusion.o ISFRunAndTumbleAndDiffusionNoLT.o ISFSwimmer.o ISFSwimmerSimp.o LaplaceTrans.o aveQ.o averSqrModTau.o clean.o ddm.o fitting.o loadFiles.o print.o shiftImage.o NILT.o NILT_RTDP.o ISFRTDP_NoLT.o ISFRTDP_NoLT_sigma.o
LFLAGS = -lm -lgsl -lgslcblas -O3 -lfftw3 -lfftw3l -fopenmp


ddm : $(Objects)
	$(CXX) -o ddm $(LFLAGS) $(Objects)

debug : $(debugObjects)
	$(CXX) -o debug $(LFLAGS) $(debugObjects)
	./debug

$(Objects) : ddm.h parameters.h
$(debugObjects) : ddm.h parameters.h
debugmain.o NILT.o NILT_RTDP.o ISFRunAndTumbleAndDiffusionNoLT.o ISFRTDP_NoLT.o : NILT.h

.PHONY : clean
clean :
	-rm $(Objects) debugmain.o