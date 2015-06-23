//
//  NILT.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "NILT.h"

NILT::NILT()
{
    CoeA=new vector<long double> [OMP_NUM_THREADS];
    for (int iter=0; iter<OMP_NUM_THREADS; ++iter)
    {
        fftwIn[iter]=fftwl_alloc_complex(M);
        fftwOut[iter]=fftwl_alloc_complex(M);
        integration[iter]=fftwl_plan_dft_1d(M, fftwIn[iter], fftwOut[iter], FFTW_FORWARD, FFTW_MEASURE);
        CoeA[iter].resize(M);
    }
}

NILT::~NILT()
{
    for (int iter=0; iter<OMP_NUM_THREADS; ++iter)
    {
        fftwl_free(fftwIn[iter]);
        fftwl_free(fftwOut[iter]);
        fftwl_destroy_plan(integration[iter]);
        CoeA[iter].clear();
    }
    delete [] CoeA;
}