//
//  NILT.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "NILT.h"
#include <omp.h>
#include <iostream>

NILT::NILT(int omp_num)
{
    OMP_NUM_THREADS=omp_num;
    integration=new fftwl_plan[omp_num];
    fftwIn=new fftwl_complex*[omp_num];
    fftwOut=new fftwl_complex*[omp_num];
    b=new long double[omp_num];
    b2=new long double[omp_num];
    sigma=new long double[omp_num];
    sigmab=new long double[omp_num];
    
    CoeA=new vector<long double> [OMP_NUM_THREADS];
    
#ifdef IfComplexIntegration
    workspace=new gsl_integration_cquad_workspace*[omp_num];
    cfun=new warper[omp_num];
    pRe=new gsl_function[omp_num];
    pIm=new gsl_function[omp_num];
#endif
    
    for (int iter=0; iter<OMP_NUM_THREADS; ++iter)
    {
        fftwIn[iter]=fftwl_alloc_complex(M);
        fftwOut[iter]=fftwl_alloc_complex(M);
        integration[iter]=fftwl_plan_dft_1d(M, fftwIn[iter], fftwOut[iter], FFTW_FORWARD, FFTW_MEASURE);
        CoeA[iter].resize(M);
        
#ifdef IfComplexIntegration
        //workspace[iter]=gsl_integration_workspace_alloc(workspaceSize);
        workspace[iter]=gsl_integration_cquad_workspace_alloc(workspaceSize);
        pRe[iter].params=&cfun[iter];
        pIm[iter].params=&cfun[iter];
        pRe[iter].function=Re;
        pIm[iter].function=Im;
#endif
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
        
#ifdef IfComplexIntegration
        //gsl_integration_workspace_free(workspace[iter]);
        gsl_integration_cquad_workspace_free(workspace[iter]);
#endif
    }
    delete [] integration;
    delete [] fftwIn;
    delete [] fftwOut;
    delete [] b;
    delete [] b2;
    delete [] sigma;
    delete [] sigmab;
    delete [] CoeA;
#ifdef IfComplexIntegration
    delete [] workspace;
    delete [] pRe;
    delete [] pIm;
    delete [] cfun;
#endif
}

void NILT::NiLT_weeks(cpx (*fun)(cpx, const long double*), const long double* para)
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
        cpx s=csigma-cb*(expterm+2.0l)/expterm;
        cpx r=-cb2*fun(s,para)/expterm;
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

double NILT::clenshaw(long double t)
{
    int tid=omp_get_thread_num();
    vector<long double>& cCoeA=CoeA[tid];
    
    long double& cb2=b2[tid];
    long double& csigmab=sigmab[tid];
    
    int idx=M-1;
    for (int iter=0; iter<M; ++iter)
    {
        if (abs(cCoeA[iter])<1e-18)
        {
            idx=iter;
            break;
        }
    }
    long double y2=0.0l;
    long double y1=cCoeA[idx];
    long double y0=0.0l;
    for (int k=idx; k>0; --k)
    {
        long double lk=(long double)k;
        y0=(2.0l*lk-1.0l-cb2*t)/lk*y1-lk/(lk+1.0l)*y2+cCoeA[k-1];
        y2=y1;
        y1=y0;
    }
    return exp(csigmab*t)*y0;
}

//Estimate parameter b by Weidemans's method, if the function is dominated by one pair of singularities \alpha1+-i\beta1. The sigma is set to be close to the singularity (sigma=alpha1+incre).
void NILT::weideman(long double alpha1, long double beta1, long double incre)
{
    long double esigma=alpha1+incre;
    long double eb=sqrt(beta1*beta1+incre*incre);
    long double eb2=eb*2;
    long double esigmab=esigma-eb;
    for (int iter=0; iter<OMP_NUM_THREADS; ++iter)
    {
        sigma[iter]=esigma;
        b[iter]=eb;
        b2[iter]=eb2;
        sigmab[iter]=esigmab;
    }
}

//Estimate parameter b by Weidemans's method, if the function is dominated by two pair of singularities \alpha1+-i\beta1. The sigma is set to be close to the singularity with largest real part, which should set to be alpha2 (sigma=alpha2+incre).
void NILT::weideman(long double alpha1, long double beta1, long double alpha2, long double beta2, long double incre)
{
    long double esigma=alpha2+incre;
    long double eb=sqrt(esigma*esigma-((esigma-alpha1)*(alpha2*alpha2+beta2*beta2)-incre*(alpha1*alpha1+beta1*beta1))/(alpha2-alpha1));
    long double eb2=eb*2;
    long double esigmab=esigma-eb;
    for (int iter=0; iter<OMP_NUM_THREADS; ++iter)
    {
        sigma[iter]=esigma;
        b[iter]=eb;
        b2[iter]=eb2;
        sigmab[iter]=esigmab;
    }
}