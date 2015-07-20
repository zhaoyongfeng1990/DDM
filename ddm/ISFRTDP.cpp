//
//  ISFRTDP_NoLT_sigma.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/7/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>

#ifdef ISFRTDP

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <vector>
#include <omp.h>

cpx ISFs(cpx s, long double* para, long double v);
cpx dvbarISFs(cpx s, long double* para, long double v);
cpx dZISFs(cpx s, long double* para, long double v);
cpx dDISFs(cpx s, long double* para, long double v);
cpx dlambdaISFs(cpx s, long double* para, long double v);

cpx ISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double vbsigma2=para[3];
    long double logfactor=para[4];
    long double vb2sigma2=para[5];
    
    long double qv=q*v;
    long double pv=exp((vb2sigma2-1)*log(v)+logfactor-v*vbsigma2);
    
    cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return atanterm/(qv-lambda*atanterm)*pv;
}

cpx dvbarISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double vbsigma2=para[3];
    long double logfactor=para[4];
    long double vb2sigma2=para[5];
    long double vb=para[6];
    long double cpsiz1=para[7];
    
    long double logv=log(v);
    long double qv=q*v;
    long double pv=exp((vb2sigma2-1)*logv+logfactor-v*vbsigma2);
    
    cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return atanterm/(qv-lambda*atanterm)*pv*vbsigma2*(1.0l-v/vb+2.0l*(logv+cpsiz1));
}

cpx dsigmaISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double vbsigma2=para[3];
    long double logfactor=para[4];
    long double vb2sigma2=para[5];
    long double vb=para[6];
    long double cpsiz1=para[7];
    long double vb2sigma3=para[8];
    
    long double logv=log(v);
    long double qv=q*v;
    long double pv=exp((vb2sigma2-1)*logv+logfactor-v*vbsigma2);
    
    cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return 2.0l*atanterm/(qv-lambda*atanterm)*pv*vb2sigma3*(v/vb-1.0l-logv-cpsiz1);
}

cpx dDISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double vbsigma2=para[3];
    long double logfactor=para[4];
    long double vb2sigma2=para[5];
    
    long double qv=q*v;
    long double pv=exp((vb2sigma2-1)*log(v)+logfactor-v*vbsigma2);
    
    long double qv2=qv*qv;
    
    cpx Dq2lambdas=Dq2lambda+s;
    cpx atanterm=atan(qv/Dq2lambdas);
    return -q*q*qv2/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm)*pv;
}

cpx dlambdaISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double vbsigma2=para[3];
    long double logfactor=para[4];
    long double vb2sigma2=para[5];
    
    long double qv=q*v;
    long double pv=exp((vb2sigma2-1)*log(v)+logfactor-v*vbsigma2);
    
    long double qv2=qv*qv;
    
    cpx Dq2lambdas=Dq2lambda+s;
    cpx atanterm=atan(qv/Dq2lambdas);
    return ((qv2+Dq2lambdas*Dq2lambdas)*atanterm*atanterm-qv2)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm)*pv;
}

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double* qArray=((dataStruct *)sdata)->q;
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    
    int num_fit=((dataStruct *)sdata)->num_fit;
    int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    int tnum_fit=num_fit*cnum_qCurve;
    
    //Get the parameters.
    long double alpha=gsl_vector_get(para, 0);
    long double vbar=gsl_vector_get(para, 1);
    long double sigma=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    
    bool breakFlag=false;
    
    double punishment=0;
    if (vbar<0)
    {
        punishment+=1e5*vbar*vbar;
        breakFlag=true;
    }
    if (sigma<0)
    {
        punishment+=1e5*sigma*sigma;
        breakFlag=true;
    }
    if (sigma>vbar)
    {
        punishment+=1e5*(sigma-vbar)*(sigma-vbar);
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
    
    //Temperary variables used for acceleration.
    long double vbsigma2=vbar/sigma/sigma;
    long double vb2sigma2=vbsigma2*vbar;
    long double logfactor=vb2sigma2*log(vbsigma2)-gsl_sf_lngamma(vb2sigma2);
    
    long double paraISF[10]={lambda, 0, 0, vbsigma2, logfactor, vb2sigma2, vbar, 0, 0, sigma};
    
    //Initialization of numerical inverse Laplace transformation solver
    int tid=omp_get_thread_num();
    
    ILT->cfun[tid].fun=ISFs;
    
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    //cout << "ok" << '\n';
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        long double A=gsl_vector_get(para, 5+2*iterqc);
        long double B=gsl_vector_get(para, 6+2*iterqc);
        
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
        
        long double q=qArray[iterqc];
        long double Dq2=D*q*q;
        double Dq2lambda=Dq2+lambda;
        
        paraISF[1]=q;
        paraISF[2]=Dq2lambda;
        
        csigma=-Dq2+0.5l;
        cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
        
        //Temperary variables used for acceleration.
        cb2=cb*2;
        csigmab=csigma-cb;
        
        //Calculate the coefficients of Laguerre polynomial series expansion.
        ILT->NiLT_weeks(paraISF);
        
        for (int iter = 0; iter<num_fit; ++iter)
        {
            int cidx=iter+iterqc*num_fit;
            long double t=tau[cidx];
            //Evaluate ISF at time t, the coefficients has been calculated.
            double rtd=ILT->clenshaw(t);
            double yi=log(A2*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B2);
            
            double result = yi - dataAry[cidx];
            
            gsl_vector_set(y, cidx, result);
        }
    }
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    int num_fit=((dataStruct *)sdata)->num_fit;
    int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    int tnum_fit=num_fit*cnum_qCurve;
    
    long double alpha=gsl_vector_get(para, 0);
    long double vbar=gsl_vector_get(para, 1);
    long double sigma=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    
    gsl_matrix_set_zero(J);
    bool breakFlag=false;
    if (vbar<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 1, 2e5 *vbar);
        }
        breakFlag=true;
    }
    if (sigma<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 2, 2e5 *sigma);
        }
        breakFlag=true;
    }
    if (sigma>vbar)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 2, 2e5 *(sigma-vbar) );
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J,iter,1) - 2e5 *(sigma-vbar) );
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
    
    double* tau=((dataStruct *)sdata)->tau;
    double* qArray=((dataStruct *)sdata)->q;
    
    //Temperary variables used for acceleration.
    long double vbsigma2=vbar/sigma/sigma;
    long double vb2sigma2=vbsigma2*vbar;
    long double logvbsigma2=log(vbsigma2);
    long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    long double vb2sigma3=vb2sigma2/sigma;
    
    long double paraISF[10]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2, vbar, cpsiz1, vb2sigma3, sigma};
    
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    NILT* dvbarILT=((dataStruct *)sdata)->dvbarISFILT;
    NILT* dsigmaILT=((dataStruct *)sdata)->dsigmaISFILT;
    NILT* dlambdaILT=((dataStruct *)sdata)->dlambdaISFILT;
    NILT* dDILT=((dataStruct *)sdata)->dDISFILT;
    
    int tid=omp_get_thread_num();
    ILT->cfun[tid].fun=ISFs;
    dvbarILT->cfun[tid].fun=dvbarISFs;
    dsigmaILT->cfun[tid].fun=dsigmaISFs;
    dlambdaILT->cfun[tid].fun=dlambdaISFs;
    dDILT->cfun[tid].fun=dDISFs;
    
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        long double A=gsl_vector_get(para, 5+2*iterqc);
        long double B=gsl_vector_get(para, 6+2*iterqc);
        
        if (A<0)
        {
            for (int iter=0; iter<tnum_fit; ++iter)
            {
                gsl_matrix_set(J, iter, 5+2*iterqc, 2e5*A);
            }
            breakFlag=true;
        }
        
        if (breakFlag)
        {
            return GSL_SUCCESS;
        }
        
        long double q=qArray[iterqc];
        long double Dq2=D*q*q;
        long double Dq2lambda=Dq2+lambda;
        
        paraISF[1]=q;
        paraISF[2]=Dq2lambda;
        
        csigma=-Dq2+0.5l;
        cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
        
        cb2=cb*2;
        csigmab=csigma-cb;
        
        dvbarILT->sigma[tid]=csigma;
        dsigmaILT->sigma[tid]=csigma;
        dlambdaILT->sigma[tid]=csigma;
        dDILT->sigma[tid]=csigma;
        
        dvbarILT->b[tid]=cb;
        dsigmaILT->b[tid]=cb;
        dlambdaILT->b[tid]=cb;
        dDILT->b[tid]=cb;
        
        dvbarILT->b2[tid]=cb2;
        dsigmaILT->b2[tid]=cb2;
        dlambdaILT->b2[tid]=cb2;
        dDILT->b2[tid]=cb2;
        
        dvbarILT->sigmab[tid]=csigmab;
        dsigmaILT->sigmab[tid]=csigmab;
        dlambdaILT->sigmab[tid]=csigmab;
        dDILT->sigmab[tid]=csigmab;
        
        //cout << "isf" << '\n';
        ILT->NiLT_weeks(paraISF);
        //cout << "dvisf" << '\n';
        dvbarILT->NiLT_weeks(paraISF);
        //cout << "dzisf" << '\n';
        dsigmaILT->NiLT_weeks(paraISF);
        //cout << "ddisf" << '\n';
        dDILT->NiLT_weeks(paraISF);
        //cout << "dlisf" << '\n';
        dlambdaILT->NiLT_weeks(paraISF);
        //cout << "ok" << '\n';
        
        for (int iter=0; iter<num_fit; ++iter)
        {
            int cidx=iter+iterqc*num_fit;
            //Temperary variables used for acceleration.
            long double t=tau[cidx];
            double rtd=ILT->clenshaw(t);
            double dvbarrtd=dvbarILT->clenshaw(t);
            double dsigmartd=dsigmaILT->clenshaw(t);
            double dlambdartd=dlambdaILT->clenshaw(t);
            double dDrtd=dDILT->clenshaw(t);
            
            double expterm=exp(-Dq2*t);
            double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
            double yi=A*dA+B;
            
            gsl_matrix_set(J, cidx, 0, A*(expterm-rtd)/yi );
            gsl_matrix_set(J, cidx, 1, -A*alpha*dvbarrtd/yi );
            gsl_matrix_set(J, cidx, 2, -A*alpha*dsigmartd/yi );
            gsl_matrix_set(J, cidx, 3, -A*alpha*dlambdartd/yi );
            gsl_matrix_set(J, cidx, 4, A*((1.0-alpha)*expterm*q*q*t-alpha*dDrtd)/yi );
            gsl_matrix_set(J, cidx, 5+2*iterqc, dA/yi );
            gsl_matrix_set(J, cidx, 6+2*iterqc, 1.0/yi );
        }
    }
    return GSL_SUCCESS;
}

#endif
