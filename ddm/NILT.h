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

const int M=512;

typedef complex<long double> cpx;

class NILT
{
public:
    NILT();
    ~NILT();
    
    void NiLT_weeks(cpx (*fun)(cpx, long double*), long double* para);
    double clenshaw(long double t);
    
    fftwl_plan integration[OMP_NUM_THREADS];
    fftwl_complex* fftwIn[OMP_NUM_THREADS];
    fftwl_complex* fftwOut[OMP_NUM_THREADS];
    vector<long double>* CoeA;
    
    long double b[OMP_NUM_THREADS];
    long double b2[OMP_NUM_THREADS];
    long double sigma[OMP_NUM_THREADS];
    long double sigmab[OMP_NUM_THREADS];
};

#endif
