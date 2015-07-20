//
//  NILT_RTDP.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/29.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "NILT.h"
#include <omp.h>
//#include <iostream>

#ifdef IfComplexIntegration

double Re(double x, void* params)
{
    warper* cfun=(warper*)params;
    //long double v=(1-x)/x;
    cpx result=cfun->fun(cfun->z, cfun->parameters, x);
    return result.real();
}

double Im(double x, void* params)
{
    warper* cfun=(warper*)params;
    
    //long double v=(1-x)/x;
    
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
    long double csigma=sigma[tid];
    long double cb=b[tid];
    long double cb2=b2[tid];
    
    for (int iter=0; iter<M; ++iter)
    {
        cpx itheta={0.0l,pi*(2.0l*iter+1.0l)/M};
        cpx expterm=exp(itheta)-1.0l;
        //cout << csigma << '\n';
        //cout << cb << '\n';
        //cout << expterm << '\n';
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
    size_t nevals;
    
    cfun[tid].parameters=para;
    cfun[tid].z=x;
    
    //gsl_integration_qagiu(&pRe[tid], 0, epsabs, epsrel, workspaceSize, workspace[tid], &realpart, &error);
    //gsl_integration_qagiu(&pIm[tid], 0, epsabs, epsrel, workspaceSize, workspace[tid], &imagpart, &error);
    
    double vb=para[6];
    double sigma=para[9];
    double lowerBound=vb-8*sigma;
    lowerBound=lowerBound<0 ? 0 : lowerBound;
    double upperBound=vb+8*sigma;
    
    gsl_integration_cquad(&pRe[tid], lowerBound, upperBound, epsabs, epsrel, workspace[tid], &realpart, &error, &nevals);
    gsl_integration_cquad(&pIm[tid], lowerBound, upperBound, epsabs, epsrel, workspace[tid], &imagpart, &error, &nevals);
    
    return cpx(realpart, imagpart);
}

#endif
