//
//  NILT.h
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef ddm_NILT_h
#define ddm_NILT_h

#include "parameters.h"
#include <fftw3.h>
#include <vector>
#include <complex>

using namespace std;

//The number of points in evaluating numerical integration in Weeks method. This also gives the number of terms in series expansion. But not all the terms in expansion is necessary.
const int M=1024;

typedef complex<long double> cpx;

#ifdef IfComplexIntegration

#include <gsl/gsl_integration.h>

//To make the iterface between real itegration and complex evaluation of functions, we use this structure to present the function that is to be integrated.
struct warper
{
    cpx z;
    cpx (*fun)(cpx z, long double* para, long double x);
    long double* parameters;
};

//Error tolerance and size of workspace for numerical integration.
const double epsabs=1e-8;
const double epsrel=1e-8;
const int workspaceSize=100000;
#endif

//Class for numerical inverse Laplace transformation.
class NILT
{
public:
    NILT(int omp_num);
    ~NILT();
    
    //Calculation of the coefficients in Laguerre polynomial expansion.
    void NiLT_weeks(cpx (*fun)(cpx, const long double*), const long double* para);
    //Clenshaw summation for function evaluation.
    double clenshaw(long double t);
    //Estimate parameter b by Weidemans's method, if the function is dominated by one pair of singularities \alpha1+-i\beta1. The sigma is set to be close to the singularity (sigma=alpha1+incre).
    void weideman(long double alpha1, long double beta1, long double incre);
    //Estimate parameter b by Weidemans's method, if the function is dominated by two pair of singularities \alpha1+-i\beta1. The sigma is set to be close to the singularity with largest real part, which should set to be alpha2 (sigma=alpha2+incre).
    void weideman(long double alpha1, long double beta1, long double alpha2, long double alpha3, long double incre);
    
    //The class is usually defined outside the parallel part of the code, to avoid allocate memory at every iteration. But this will causs memory conflict in shared memory model like openMP. So everything should be kept as a list with number of elements equals to number of threads, and different thread uses different element.
    
    int OMP_NUM_THREADS;    //Number of threads
    
    //FFTW stuff used in calculating numerical integration.
    fftwl_plan* integration;    //Integration is done by FFT.
    fftwl_complex** fftwIn;
    fftwl_complex** fftwOut;
    //Coefficients in Laguerre polynomial expansion.
    vector<long double>* CoeA;
    
    //Important paramters in solver.
    long double* b;
    long double* sigma;
    
    //b2 = b*2
    long double* b2;
    //sigmab = sigma-b
    long double* sigmab;
    
#ifdef IfComplexIntegration
    //Calculation of the coefficients in Laguerre polynomial expansion, with integration on v.
    void NiLT_weeks(long double* para);
    //Numerical evaluation of function to be inverse transformed.
    cpx invfun(cpx x, long double* para);
    
    //Pointers of Re and Im
    gsl_function* pRe;
    gsl_function* pIm;
    
    //gsl_integration_workspace* workspace;
    gsl_integration_cquad_workspace** workspace;
    
    //To make the iterface between real itegration and complex evaluation of functions, we use this structure to present the function that is to be integrated.
    warper* cfun;
#endif
};

#ifdef IfComplexIntegration
//Real part and image part of the complex function. Defined outside the class to make use of their pointers.
double Re(double x, void* params);
double Im(double x, void* params);
#endif

#endif
