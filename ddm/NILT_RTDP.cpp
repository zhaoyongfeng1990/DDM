//
//  NILT_RTDP.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/29.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "NILT.h"
#include <omp.h>

#ifdef IfComplexIntegration

double Re(double x, void* params)
{
    warper* cfun=(warper*)params;
    cpx result=cfun->fun(cfun->z, cfun->parameters, x);
    return result.real();
}

double Im(double x, void* params)
{
    warper* cfun=(warper*)params;
    cpx result=cfun->fun(cfun->z, cfun->parameters, x);
    return result.imag();
}

void NILT::NiLT_weeks(long double* para)
{
    int tid=omp_get_thread_num();
    fftwl_complex* cfftwIn=fftwIn[tid];
    fftwl_complex* cfftwOut=fftwOut[tid];
    fftwl_plan& cPlan=integration[tid];
    vector<long double>& cCoeA=CoeA[tid];
    long double& csigma=sigma[tid];
    long double& cb=b[tid];
    long double& cb2=b2[tid];
    
    for (int iter=0; iter<M; ++iter)
    {
        cpx itheta={0.0l,pi*(2.0l*iter+1.0l)/M};
        cpx expterm=exp(itheta)-1.0l;
        cpx s=csigma-cb*(expterm+2.0l)/expterm;
        cpx r=-cb2*invfun(s,para)/expterm;
        cfftwIn[iter][0]=real(r);
        cfftwIn[iter][1]=imag(r);
    }
    fftwl_execute(cPlan);
    for (int iter=0; iter<M; ++iter)
    {
        cpx temp={cfftwOut[iter][0], cfftwOut[iter][1]};
        cpx itheta={0,-iter*pi/M};
        temp*=exp(itheta);
        cCoeA[iter]=real(temp)/M;
    }
}

cpx NILT::invfun(cpx x, long double* para)
{
    int tid=omp_get_thread_num();
    double realpart, imagpart, error;
    
    cfun[tid].parameters=para;
    cfun[tid].z=x;
    
    gsl_integration_qagiu(&pRe[tid], 0, epsabs, epsrel, workspaceSize, workspace[tid], &realpart, &error);
    gsl_integration_qagiu(&pIm[tid], 0, epsabs, epsrel, workspaceSize, workspace[tid], &imagpart, &error);
    
    return cpx(realpart, imagpart);
}

#endif
