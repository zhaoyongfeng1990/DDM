//
//  NILT_RTDP.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/29.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "NILT.h"
#include <omp.h>
#include <iostream>

#ifdef IfComplexIntegration

//Real part and image part of the complex function. Defined outside the class to make use of their pointers.
double Re(double x, void* params)
{
    warper* cfun=(warper*)params;
    long double v=(1-x)/x;  //Transform the integration range from [0,Infinity] to [0,1]
    //cpx result=cfun->fun(cfun->z, cfun->parameters, x);
    cpx result=cfun->fun(cfun->z, cfun->parameters, v);
    return result.real()/x/x;
}

double Im(double x, void* params)
{
    warper* cfun=(warper*)params;
    long double v=(1-x)/x;  //Transform the integration range from [0,Infinity] to [0,1]
    //cpx result=cfun->fun(cfun->z, cfun->parameters, x);
    cpx result=cfun->fun(cfun->z, cfun->parameters, v);
    return result.imag()/x/x;
}

//Calculation of the coefficients in Laguerre polynomial expansion, with integration on v.
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

//Numerical evaluation of function to be inverse transformed.
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
    double lowerBound=0;//vb-8*sigma;
    //lowerBound=lowerBound<0 ? 0 : lowerBound;
    double upperBound=1;//vb+8*sigma;
    
    gsl_integration_cquad(&pRe[tid], lowerBound, upperBound, epsabs, epsrel, workspace[tid], &realpart, &error, &nevals);
    gsl_integration_cquad(&pIm[tid], lowerBound, upperBound, epsabs, epsrel, workspace[tid], &imagpart, &error, &nevals);
    
    return cpx(realpart, imagpart);
}

int NILT::optimize_incre(long double alpha1, long double beta1, long double* para, long double beginTime, long double finalTime)
{
    const long double golden=(sqrt(5.0l)-1.0l)/2.0l;
    long double left=0.0l;
    long double right=abs(2.0l*alpha1);
    long double mid=left+golden*(right-left);
    bool ifMidRight=1;
    bool ifNeedRecalculate=0;
    //weideman(alpha1, beta1, left);
    //NiLT_weeks(para);
    //long double left_value=estimate_Err(beginTime, finalTime);
    weideman(alpha1, beta1, mid);
    NiLT_weeks(para);
    long double mid_value=estimate_Err(beginTime, finalTime);
    //weideman(alpha1, beta1, right);
    //NiLT_weeks(para);
    //long double right_value=estimate_Err(beginTime, finalTime);
    while (right-left>1e-5)
    {
        if (ifMidRight)
        {
            long double test=left+(mid-left)*golden;
            weideman(alpha1, beta1, test);
            NiLT_weeks(para);
            long double test_value=estimate_Err(beginTime, finalTime);
            if (test_value>mid_value)
            {
                left=test;
                //left_value=test_value;
                ifMidRight=0;
                ifNeedRecalculate=1;
            }
            else
            {
                right=mid;
                //right_value=mid_value;
                mid=test;
                mid_value=test_value;
            }
        }
        else
        {
            long double test=right-(right-mid)*golden;
            weideman(alpha1, beta1, test);
            NiLT_weeks(para);
            long double test_value=estimate_Err(beginTime, finalTime);
            if (test_value>mid_value)
            {
                right=test;
                //right_value=test_value;
                ifMidRight=1;
                ifNeedRecalculate=1;
            }
            else
            {
                left=mid;
                //left_value=mid_value;
                mid=test;
                mid_value=test_value;
            }
        }
    }
    if (ifNeedRecalculate)
    {
        weideman(alpha1, beta1, mid);
        NiLT_weeks(para);
    }
    if(estimate_Err(beginTime, finalTime)>tol_NiLT)
        return 0;
    return 1;
}

int NILT::optimize_incre(long double alpha1, long double beta1, long double alpha2, long double beta2, long double* para, long double beginTime, long double finalTime)
{
    const long double golden=(sqrt(5.0l)-1.0l)/2.0l;
    long double left=0.0l;
    long double right=abs(2.0l*alpha2);
    long double mid=left+golden*(right-left);
    bool ifMidRight=1;
    bool ifNeedRecalculate=0;
    //weideman(alpha1, beta1, alpha2, beta2, left);
    //NiLT_weeks(para);
    //long double left_value=estimate_Err(beginTime, finalTime);
    weideman(alpha1, beta1, alpha2, beta2, mid);
    NiLT_weeks(para);
    long double mid_value=estimate_Err(beginTime, finalTime);
    //weideman(alpha1, beta1, alpha2, beta2, right);
    //NiLT_weeks(para);
    //long double right_value=estimate_Err(beginTime, finalTime);
    while (right-left>1e-5)
    {
        //cout << mid << endl;
        if (ifMidRight)
        {
            long double test=left+(mid-left)*golden;
            weideman(alpha1, beta1, alpha2, beta2, test);
            NiLT_weeks(para);
            long double test_value=estimate_Err(beginTime, finalTime);
            if (test_value>mid_value)
            {
                left=test;
                //left_value=test_value;
                ifMidRight=0;
                ifNeedRecalculate=1;
            }
            else
            {
                right=mid;
                //right_value=mid_value;
                mid=test;
                mid_value=test_value;
            }
        }
        else
        {
            long double test=right-(right-mid)*golden;
            weideman(alpha1, beta1, alpha2, beta2, test);
            NiLT_weeks(para);
            long double test_value=estimate_Err(beginTime, finalTime);
            if (test_value>mid_value)
            {
                right=test;
                //right_value=test_value;
                ifMidRight=1;
                ifNeedRecalculate=1;
            }
            else
            {
                left=mid;
                //left_value=mid_value;
                mid=test;
                mid_value=test_value;
            }
        }
    }
    if (ifNeedRecalculate)
    {
        weideman(alpha1, beta1, alpha2, beta2, mid);
        NiLT_weeks(para);
    }
    if(estimate_Err(beginTime, finalTime)>tol_NiLT)
        return 0;
    return 1;
}

#endif
