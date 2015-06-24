//
//  ISFRunAndTumbleAndDiffusionNoLT.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/18.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"

#ifdef ISFRunAndTumbleAndDiffusionNoLT

#include <vector>
#include <omp.h>

cpx ISFs(cpx s, long double* para);
cpx dvISFs(cpx s, long double* para);
cpx dDISFs(cpx s, long double* para);
cpx dlambdaISFs(cpx s, long double* para);

cpx ISFs(cpx s, long double* para)
{
    long double qv=para[0];
    long double lambda=para[1];
    long double Dq2=para[2];
    long double q=para[3];
    long double Dq2lambda=para[4];
    cpx atanterm=atan(qv/(Dq2lambda+s));
    return atanterm/(qv-lambda*atanterm);
}

cpx dvISFs(cpx s, long double* para)
{
    long double qv=para[0];
    long double lambda=para[1];
    long double Dq2=para[2];
    long double q=para[3];
    long double Dq2lambda=para[4];
    long double qv2=para[5];
    cpx Dq2lambdas=Dq2lambda+s;
    cpx atanterm=atan(qv/Dq2lambdas);
    return q*(qv*Dq2lambdas-(qv2+Dq2lambdas*Dq2lambdas)*atanterm)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm);
}

cpx dDISFs(cpx s, long double* para)
{
    long double qv=para[0];
    long double lambda=para[1];
    long double Dq2=para[2];
    long double q=para[3];
    long double Dq2lambda=para[4];
    long double qv2=para[5];
    cpx Dq2lambdas=Dq2lambda+s;
    cpx atanterm=atan(qv/Dq2lambdas);
    return -q*q*qv2/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm);
}

cpx dlambdaISFs(cpx s, long double* para)
{
    long double qv=para[0];
    long double lambda=para[1];
    long double Dq2=para[2];
    long double q=para[3];
    long double Dq2lambda=para[4];
    long double qv2=para[5];
    cpx Dq2lambdas=Dq2lambda+s;
    cpx atanterm=atan(qv/Dq2lambdas);
    return ((qv2+Dq2lambdas*Dq2lambdas)*atanterm*atanterm-qv2)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm);
}

void NILT::NiLT_weeks(cpx (*fun)(cpx, long double*), long double* para)
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
        if (abs(cCoeA[iter])<1e-16)
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

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    long double q=((dataStruct *)sdata)->q;
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    
    //Get the parameters.
    long double alpha=gsl_vector_get(para, 0);
    long double v0=gsl_vector_get(para, 1);
    long double lambda=gsl_vector_get(para, 2);
    long double D=gsl_vector_get(para, 3);
    long double A=gsl_vector_get(para, 4);
    long double B=gsl_vector_get(para, 5);
    
    //cout << lambda << endl;

    long double kv0=q*v0;
    long double Dq2=D*q*q;
    long double Dq2lambda=Dq2+lambda;
    long double paraISF[5]={kv0, lambda, Dq2, q, Dq2lambda};
    
    long double qvlambda=kv0/lambda;
    
    int tid=omp_get_thread_num();
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    if (qvlambda>(pi/2))
    {
        csigma=-Dq2-lambda+0.01l;
        long double temp=(csigma+Dq2lambda);
        cb=sqrt(kv0*kv0-temp*temp);
    }
    else
    {
        long double alpha21=kv0/tan(qvlambda);
        csigma=alpha21-Dq2lambda+0.01l;
        long double alpha1=-Dq2lambda;
        long double alpha2=alpha21+alpha1;
        cb=sqrt(csigma*csigma-((csigma-alpha1)*alpha2*alpha2-(csigma-alpha2)*(alpha1*alpha1+kv0*kv0))/alpha21);
    }
    cb2=cb*2;
    csigmab=csigma-cb;
    
    ILT->NiLT_weeks(ISFs, paraISF);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        long double t=tau[iter];
        double rtd=ILT->clenshaw(t);
        //Temperary variables used for acceleration.
        double yi=(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
        
        double result = yi - dataAry[iter];
        //Punishment terms, to make constrains in parameter space.
        if (v0<0)
        {
            result += 1e5*v0*v0;
        }
        if (lambda<0)
        {
            result += 1e5*lambda*lambda;
        }
        if (D<0)
        {
            result += 1e5*D*D;
        }
        if (alpha<0)
        {
            result += 1e5*alpha*alpha;
        }
        if (alpha>1)
        {
            result += 1e5*(alpha-1)*(alpha-1);
        }
        if (A<0)
        {
            result += 1e5*A*A;
        }
        
        gsl_vector_set(y, iter, result);
    }
    
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    long double q=((dataStruct *)sdata)->q;
    
    long double alpha=gsl_vector_get(para, 0);
    long double v0=gsl_vector_get(para, 1);
    long double lambda=gsl_vector_get(para, 2);
    long double D=gsl_vector_get(para, 3);
    long double A=gsl_vector_get(para, 4);
    long double B=gsl_vector_get(para, 5);
    
    long double kv0=q*v0;
    long double Dq2=D*q*q;
    long double Dq2lambda=Dq2+lambda;
    long double qv2=kv0*kv0;
    long double paraISF[6]={kv0, lambda, Dq2, q, Dq2lambda, qv2};
    
    long double qvlambda=q*v0/lambda;
    
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    NILT* dvILT=((dataStruct *)sdata)->dvISFILT;
    NILT* dlambdaILT=((dataStruct *)sdata)->dlambdaISFILT;
    NILT* dDILT=((dataStruct *)sdata)->dDISFILT;
    
    int tid=omp_get_thread_num();
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    if (qvlambda>(pi/2))
    {
        csigma=-Dq2-lambda+0.01l;
        long double temp=(csigma+Dq2lambda);
        cb=sqrt(kv0*kv0-temp*temp);
    }
    else
    {
        long double alpha21=kv0/tan(qvlambda);
        csigma=alpha21-Dq2lambda+0.01l;
        long double alpha1=-Dq2lambda;
        long double alpha2=alpha21+alpha1;
        cb=sqrt(csigma*csigma-((csigma-alpha1)*alpha2*alpha2-(csigma-alpha2)*(alpha1*alpha1+kv0*kv0))/alpha21);
    }
    cb2=cb*2;
    csigmab=csigma-cb;
    
    dvILT->sigma[tid]=csigma;
    dlambdaILT->sigma[tid]=csigma;
    dDILT->sigma[tid]=csigma;
    
    dvILT->b[tid]=cb;
    dlambdaILT->b[tid]=cb;
    dDILT->b[tid]=cb;
    
    dvILT->b2[tid]=cb2;
    dlambdaILT->b2[tid]=cb2;
    dDILT->b2[tid]=cb2;
    
    dvILT->sigmab[tid]=csigmab;
    dlambdaILT->sigmab[tid]=csigmab;
    dDILT->sigmab[tid]=csigmab;
    
    ILT->NiLT_weeks(ISFs, paraISF);
    dvILT->NiLT_weeks(dvISFs, paraISF);
    dDILT->NiLT_weeks(dDISFs, paraISF);
    dlambdaILT->NiLT_weeks(dlambdaISFs, paraISF);
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        long double t=tau[iter];
        double rtd=ILT->clenshaw(t);
        double dvrtd=dvILT->clenshaw(t);
        double dlambdartd=dlambdaILT->clenshaw(t);
        double dDrtd=dDILT->clenshaw(t);
        
        double dA=(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd);
        double yi=A*dA+B;
        
        gsl_matrix_set(J, iter, 0, A*(exp(-Dq2*t)-rtd) );
        gsl_matrix_set(J, iter, 1, -A*alpha*dvrtd );
        
        gsl_matrix_set(J, iter, 2, -A*alpha*dlambdartd );
        
        gsl_matrix_set(J, iter, 3, A*((alpha-1.0)*exp(-Dq2*t)*q*q*t-alpha*dDrtd) );
        
        gsl_matrix_set(J, iter, 4, dA);
        gsl_matrix_set(J, iter, 5, 1.0 );
        
        //Punishment terms
        if (alpha<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *alpha);
        }
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
        }
        if (v0<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *v0);
        }
        if (lambda<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *lambda);
        }
        if (D<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *D);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif
