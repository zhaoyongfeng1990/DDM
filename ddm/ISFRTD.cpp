//
//  ISFRTD.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/18.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"
#include <iostream>

#ifdef ISFRTD

#include <vector>
#include <omp.h>

cpx ISFs(cpx s, const long double* para);
cpx dvISFs(cpx s, const long double* para);
cpx dDISFs(cpx s, const long double* para);
cpx dlambdaISFs(cpx s, const long double* para);

cpx ISFs(cpx s, const long double* para)
{
    //Temperary variables used for acceleration.
    const long double qv=para[0];
    const long double lambda=para[1];
    const long double Dq2=para[2];
    const long double q=para[3];
    const long double Dq2lambda=para[4];
    const cpx atanterm=atan(qv/(Dq2lambda+s));
    return atanterm/(qv-lambda*atanterm);
}

cpx dvISFs(cpx s, const long double* para)
{
    //Temperary variables used for acceleration.
    const long double qv=para[0];
    const long double lambda=para[1];
    const long double Dq2=para[2];
    const long double q=para[3];
    const long double Dq2lambda=para[4];
    const long double qv2=para[5];
    const cpx Dq2lambdas=Dq2lambda+s;
    const cpx atanterm=atan(qv/Dq2lambdas);
    return q*(qv*Dq2lambdas-(qv2+Dq2lambdas*Dq2lambdas)*atanterm)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm);
}

cpx dDISFs(cpx s, const long double* para)
{
    //Temperary variables used for acceleration.
    const long double qv=para[0];
    const long double lambda=para[1];
    const long double Dq2=para[2];
    const long double q=para[3];
    const long double Dq2lambda=para[4];
    const long double qv2=para[5];
    const cpx Dq2lambdas=Dq2lambda+s;
    const cpx atanterm=atan(qv/Dq2lambdas);
    return -q*q*qv2/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm);
}

cpx dlambdaISFs(cpx s, const long double* para)
{
    //Temperary variables used for acceleration.
    const long double qv=para[0];
    const long double lambda=para[1];
    const long double Dq2=para[2];
    const long double q=para[3];
    const long double Dq2lambda=para[4];
    const long double qv2=para[5];
    const cpx Dq2lambdas=Dq2lambda+s;
    const cpx atanterm=atan(qv/Dq2lambdas);
    return ((qv2+Dq2lambdas*Dq2lambdas)*atanterm*atanterm-qv2)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm);
}

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    const double* dataAry=((dataStruct *)sdata)->data;
    const double* tau=((dataStruct *)sdata)->tau;
    const double* qArray=((dataStruct *)sdata)->q;
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    
    const int num_fit=((dataStruct *)sdata)->num_fit;
    const int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    const int tnum_fit=num_fit*cnum_qCurve;
    
    //Get the parameters.
    const long double alpha=gsl_vector_get(para, 0);
    const long double v0=gsl_vector_get(para, 1);
    const long double lambda=gsl_vector_get(para, 2);
    const long double D=gsl_vector_get(para, 3);
    
    bool breakFlag=false;
    
    double punishment=0;
    if (v0<0)
    {
        punishment+=1e5*v0*v0;
        breakFlag=true;
    }
    if (lambda<0)
    {
        punishment+=1e5*lambda*lambda;
        breakFlag=true;
    }
    if (D<0)
    {
        punishment+=1e5*D*D;
        breakFlag=true;
    }
    if (D>10)
    {
        punishment+=1e5*(D-10)*(D-10);
        breakFlag=true;
    }
    if (alpha<0)
    {
        punishment+=1e5*alpha*alpha;
        breakFlag=true;
    }
    if (alpha>1)
    {
        punishment+=1e5*(alpha-1)*(alpha-1);
        breakFlag=true;
    }
    
    if (breakFlag)
    {
        for (int iter = 0; iter<tnum_fit; ++iter)
        {
            gsl_vector_set(y, iter, punishment - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    
    //Initialization of numerical inverse Laplace transformation solver
    const int tid=omp_get_thread_num();
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    //Temperary variables used for acceleration.
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        const long double A=gsl_vector_get(para, 4+2*iterqc);
        const long double B=gsl_vector_get(para, 5+2*iterqc);
        
        if (A<0)
        {
            punishment+=1e5*A*A;
            breakFlag=true;
        }
        if (breakFlag)
        {
            for (int iter = 0; iter<tnum_fit; ++iter)
            {
                gsl_vector_set(y, iter, punishment - dataAry[iter]);
            }
            return GSL_SUCCESS;
        }
        
        const long double q=qArray[iterqc];
        const long double kv0=q*v0;
        const long double Dq2=D*q*q;
        const long double Dq2lambda=Dq2+lambda;
        const long double qvlambda=kv0/lambda;
        const long double paraISF[5]={kv0, lambda, Dq2, q, Dq2lambda};
        //Set sigma and b in iLT solver
        if (qvlambda>(pi/2))
        {
            csigma=-Dq2lambda+1.0l;
            cb=sqrt(kv0*kv0+1.0l*1.0l);
        }
        else
        {
            const long double alpha21=kv0/tan(qvlambda);
            const long double alpha1=-Dq2lambda;
            const long double alpha2=alpha21+alpha1;
            csigma=alpha2+0.1l;
            cb=sqrt(csigma*csigma-((csigma-alpha1)*alpha2*alpha2-0.1l*(alpha1*alpha1+kv0*kv0))/alpha21);
        }
        //Temperary variables used for acceleration.
        cb2=cb*2;
        csigmab=csigma-cb;
        
        //Calculate the coefficients of Laguerre polynomial series expansion.
        ILT->NiLT_weeks(ISFs, paraISF);
        
        for (int iter = 0; iter<num_fit; ++iter)
        {
            const int cidx=iter+iterqc*num_fit;
            const long double t=tau[cidx];
            //Evaluate ISF at time t, the coefficients has been calculated.
            const double rtd=ILT->clenshaw(t);
            const double yi=(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
            const double weight=1.0/sqrt(dataAry[cidx]);
            
            double result = (yi - dataAry[cidx])*weight;
            
            gsl_vector_set(y, cidx, result);
        }
    }
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    const int num_fit=((dataStruct *)sdata)->num_fit;
    const int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    const int tnum_fit=num_fit*cnum_qCurve;
    
    const long double alpha=gsl_vector_get(para, 0);
    const long double v0=gsl_vector_get(para, 1);
    const long double lambda=gsl_vector_get(para, 2);
    const long double D=gsl_vector_get(para, 3);
    
    gsl_matrix_set_zero(J);
    bool breakFlag=false;
    if (v0<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 1, 2e5 *v0);
        }
        breakFlag=true;
    }
    if (lambda<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 3, 2e5*lambda);
        }
        breakFlag=true;
    }
    if (alpha<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 2e5*alpha);
        }
        breakFlag=true;
    }
    if (alpha>1)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 2e5*(alpha-1));
        }
        breakFlag=true;
    }
    if (D<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 4, 2e5*D);
        }
        breakFlag=true;
    }
    if (D>10)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 4, 2e5*(D-10));
        }
        breakFlag=true;
    }
    
    if (breakFlag)
    {
        return GSL_SUCCESS;
    }
    
    const double* tau=((dataStruct *)sdata)->tau;
    const double* qArray=((dataStruct *)sdata)->q;
    const double* dataAry=((dataStruct *)sdata)->data;
    
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    NILT* dvILT=((dataStruct *)sdata)->dvISFILT;
    NILT* dlambdaILT=((dataStruct *)sdata)->dlambdaISFILT;
    NILT* dDILT=((dataStruct *)sdata)->dDISFILT;
    
    const int tid=omp_get_thread_num();
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        const long double A=gsl_vector_get(para, 4+2*iterqc);
        const long double B=gsl_vector_get(para, 5+2*iterqc);
        
        if (A<0)
        {
            for (int iter=0; iter<tnum_fit; ++iter)
            {
                gsl_matrix_set(J, iter, 4+2*iterqc, 2e5*A);
            }
            breakFlag=true;
        }
        
        if (breakFlag)
        {
            return GSL_SUCCESS;
        }
        
        const long double q=qArray[iterqc];
        //Temperary variables used for acceleration.
        const long double kv0=q*v0;
        const long double Dq2=D*q*q;
        const long double Dq2lambda=Dq2+lambda;
        const long double qv2=kv0*kv0;
        const long double paraISF[6]={kv0, lambda, Dq2, q, Dq2lambda, qv2};
        
        const long double qvlambda=q*v0/lambda;
        
        if (qvlambda>(pi/2))
        {
            csigma=-Dq2lambda+1.0l;
            cb=sqrt(kv0*kv0+1.0l*1.0l);
        }
        else
        {
            const long double alpha21=kv0/tan(qvlambda);
            const long double alpha1=-Dq2lambda;
            const long double alpha2=alpha21+alpha1;
            csigma=alpha2+0.1l;
            cb=sqrt(csigma*csigma-((csigma-alpha1)*alpha2*alpha2-0.1l*(alpha1*alpha1+kv0*kv0))/alpha21);
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
            const int cidx=iter+iterqc*num_fit;
            //Temperary variables used for acceleration.
            const long double t=tau[cidx];
            const double rtd=ILT->clenshaw(t);
            const double dvrtd=dvILT->clenshaw(t);
            const double dlambdartd=dlambdaILT->clenshaw(t);
            const double dDrtd=dDILT->clenshaw(t);
            
            const double expterm=exp(-Dq2*t);
            const double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
            //double yi=A*dA+B;
            const double weight=1.0/sqrt(dataAry[cidx]);
            
            gsl_matrix_set(J, cidx, 0, A*(expterm-rtd)*weight );
            gsl_matrix_set(J, cidx, 1, -A*alpha*dvrtd*weight );
            
            gsl_matrix_set(J, cidx, 2, -A*alpha*dlambdartd*weight );
            
            gsl_matrix_set(J, cidx, 3, A*((1.0-alpha)*expterm*q*q*t-alpha*dDrtd)*weight );
            
            gsl_matrix_set(J, cidx, 4, dA*weight );
            gsl_matrix_set(J, cidx, 5, 1.0*weight );
        }
    }
    return GSL_SUCCESS;
}
#endif
