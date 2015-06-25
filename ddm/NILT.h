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

const int M=1024;
//The number of points in evaluating numerical integration in Weeks method. This also gives the number of terms in series expansion. But not all the terms in expansion is necessary.

typedef complex<long double> cpx;

//Class for numerical inverse Laplace transformation.
class NILT
{
public:
    NILT();
    ~NILT();
    
    void NiLT_weeks(cpx (*fun)(cpx, long double*), long double* para);
    //Calculation of the coefficients in Laguerre polynomial expansion.
    double clenshaw(long double t);
    //Clenshaw summation for function evaluation.
    void weideman(long double alpha1, long double beta1, long double incre);
    //Estimate parameter b by Weidemans's method, if the function is dominated by one pair of singularities \alpha1+-i\beta1. The sigma is set to be close to the singularity (sigma=alpha1+incre).
    void weideman(long double alpha1, long double beta1, long double alpha2, long double alpha3, long double incre);
    //Estimate parameter b by Weidemans's method, if the function is dominated by two pair of singularities \alpha1+-i\beta1. The sigma is set to be close to the singularity with largest real part, which should set to be alpha2 (sigma=alpha2+incre).
    
    
    fftwl_plan integration[OMP_NUM_THREADS];
    fftwl_complex* fftwIn[OMP_NUM_THREADS];
    fftwl_complex* fftwOut[OMP_NUM_THREADS];
    //FFTW stuff used in calculating numerical integration.
    vector<long double>* CoeA;
    //Coefficients in Laguerre polynomial expansion.
    
    long double b[OMP_NUM_THREADS];
    long double sigma[OMP_NUM_THREADS];
    //Important paramters in solver.
    
    long double b2[OMP_NUM_THREADS];
    //b2 is b*2
    long double sigmab[OMP_NUM_THREADS];
    //sigmab is sigma-b
};

#endif
