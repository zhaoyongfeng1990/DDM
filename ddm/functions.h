//
//  functions.h
//  ddm
//
//  Created by Zhao Yongfeng on 14/11/3.
//  Copyright (c) 2014年 ZYF. All rights reserved.
//

#ifndef __ddm__functions__
#define __ddm__functions__

//If the platform is Windows, uncomment this line.
//#define WINDOWS
#ifdef WINDOWS
#pragma comment(lib,"libgsl.lib")
#pragma comment(lib, "libgslcblas.lib")
#pragma comment(lib, "libfftw3-3.lib")
#endif

///////////////////

//Macro for speeding up GSL
#define HAVE_INLINE

//Switch of different models
//#define ISFRUNANDTUMBLE
//#define ISFSWIMMER
//#define ISFRUNANDTUMBLE_3D
#define ISFRunAndTumbleAndDiffusion

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <fftw3.h>	//FFTW
#include <fstream>
#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <omp.h>	//OpenMP support

#ifdef WINDOWS
#include <cstdint>
#endif

using namespace std;

//The data structures for reading tiff
typedef struct tagTIFFILEHEAD{
    unsigned short biOrder;
    unsigned short version;
    uint32_t offsetIFD;
} TIFFILEHEAD;

typedef struct tagTIFFDE{
    unsigned short tag;
    unsigned short type;
    uint32_t length;
    uint32_t valueOffset;
} TIFFDE;

typedef struct tagTIFFIFD{
    unsigned short DEC;
    TIFFDE* DE;
    uint32_t offset;
    tagTIFFIFD* nextIFD;
} TIFFIFD;

gsl_matrix* readTiff(const string tifName);     //Reading tiff images.
gsl_matrix* readSim(const string fileName);     //Reading simulation txt files.
void shiftImage(gsl_matrix** img, const int drow, const int dcol);      //Used for aligning viberating images.
double shiftCorrelation2d(const gsl_matrix* img1, const gsl_matrix* img2, const int drow, const int dcol, const gsl_vector* aveVec);
//Calculating shifted correlation for two images, drow and dcol indicate the displacement between two images.
double correlation2d(const gsl_matrix* img1, const gsl_matrix* img2, const gsl_vector* aveVec);
//Calculating correlation for two images.
int find(const vector<int>& vec, const int value);
//Find the position of a particular member in a vector.
int quickFind(const vector<int>& vec, const int value);
//Find the position of a particular member in a sorted vector with increasing order, using dichotomy.

const int dim = 1024;    //Size of the image
const int numOfSeq = 4500;  //Number of total time points in experiment
const int numOfDiff = 4000; //Number of tau
const double dx = 6.5 / 4.0; // 0.65;   //Pixel size
const double pi = 3.14159265358979323846264338327950288419716939937510;
const int winDim = 5;   //Size of searching window. Used in aligning images.
const int dimk = dim / 2 + 1;   //Number of wavenumber
const double dq = 2 * pi / dx / dim;    //q step after FFT
const double qmax=sqrt(dim*dim/2)*dq;   //Maximum possible value of q
const double dt = 1.0 / 100;    //Time step
const int numOfk = dim*dimk;    //Number of data points after FFT
const int num_fit = numOfDiff;  //Number of data points used in fitting.
const double sqrtpi = sqrt(pi); //Constant for convenience
const double precision = 1e-15; //Used in 2D R&T model. The precision of numerical evaluation.
const int maxIter=20000;    //Maximum iteration number in fitting

const double qstep=0.01;    //The width of cirque when averaging the direction of q.

//Data stuct used in GSL fitting algorithm.
typedef struct
{
	double* data;
	double* tau;
	double q;
} dataStruct;

//Number of parameters
#ifdef ISFSWIMMER
const int numOfPara = 6;
#endif
#ifdef ISFRUNANDTUMBLE
const int numOfPara = 4;
#endif
#ifdef ISFRUNANDTUMBLE_3D
const int numOfPara=4;
#define NeedLaplaceTrans
#endif
#ifdef ISFRunAndTumbleAndDiffusion
const int numOfPara=6; 
#define NeedLaplaceTrans
#endif


//ISF and its Jacobian. Used in GSL fitting algorithm.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y);
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J);
int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J);
//Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.
int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole);
//Test fitting result using error estimation from covariance matrix, not reliable. tol is the relative tolerance of error.
int covar_rel_test(const gsl_matrix* J, const gsl_vector* x, double tol);


#endif /* defined(__ddm__functions__) */
