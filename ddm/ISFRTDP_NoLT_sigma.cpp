//
//  ISFRTDP_NoLT_sigma.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/7/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>

#ifdef ISFRTDPNoLT_sigma

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

#ifndef MultiQFit
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    long double q=((dataStruct *)sdata)->q;
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    
    int num_fit=((dataStruct *)sdata)->num_fit;
    
    //Get the parameters.
    long double alpha=gsl_vector_get(para, 0);
    long double vbar=gsl_vector_get(para, 1);
    long double sigma=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    long double A=gsl_vector_get(para, 5);
    long double B=gsl_vector_get(para, 6);
    
    long double Dq2=D*q*q;
    
    if (sigma<0)
    {
        for (int iter = 0; iter<num_fit; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*sigma*sigma - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    if (vbar<0)
    {
        for (int iter = 0; iter<num_fit; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*vbar*vbar - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    if (lambda<0)
    {
        for (int iter = 0; iter<num_fit; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*lambda*lambda - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }

    
    ////cout << lambda << endl;
    
    //Temperary variables used for acceleration.
    long double Dq2lambda=Dq2+lambda;
    long double vbsigma2=vbar/sigma/sigma;
    long double vb2sigma2=vbsigma2*vbar;
    long double logfactor=vb2sigma2*log(vbsigma2)-gsl_sf_lngamma(vb2sigma2);
    
    long double paraISF[6]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2};
    
    //Initialization of numerical inverse Laplace transformation solver
    int tid=omp_get_thread_num();
    
//    if (tid==0)
//    {
//        cout << "debug debug \n";
//        cout << q << endl;
//        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
//        {
//            cout << gsl_vector_get(para, iterpara) << '\n';
//        }
//        cout << '\n';
//        cout << "debug debug \n";
//    }
    
    ILT->cfun[tid].fun=ISFs;
    
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    //Set sigma and b in iLT solver
    csigma=-Dq2+0.5l;
    cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
    
    //Temperary variables used for acceleration.
    cb2=cb*2;
    csigmab=csigma-cb;
    
    //Calculate the coefficients of Laguerre polynomial series expansion.
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "ok" << endl;
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        long double t=tau[iter];
        //Evaluate ISF at time t, the coefficients has been calculated.
        double rtd=ILT->clenshaw(t);
        double yi=log(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
        
        double result = yi - dataAry[iter];
        //Punishment terms, to make constrains in parameter space.
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
    
    int num_fit=((dataStruct *)sdata)->num_fit;
    
    long double alpha=gsl_vector_get(para, 0);
    long double vbar=gsl_vector_get(para, 1);
    long double sigma=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    long double A=gsl_vector_get(para, 5);
    long double B=gsl_vector_get(para, 6);
    
    long double Dq2=D*q*q;
    
    if (sigma<0)
    {
        for (int iter=0; iter<num_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 0);
            gsl_matrix_set(J, iter, 2, 2e5 *sigma);
            gsl_matrix_set(J, iter, 3, 0);
            gsl_matrix_set(J, iter, 4, 0);
            gsl_matrix_set(J, iter, 5, 0);
            gsl_matrix_set(J, iter, 6, 0);
        }
        return GSL_SUCCESS;
    }
    if (vbar<0)
    {
        for (int iter=0; iter<num_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 2e5 *vbar);
            gsl_matrix_set(J, iter, 2, 0);
            gsl_matrix_set(J, iter, 3, 0);
            gsl_matrix_set(J, iter, 4, 0);
            gsl_matrix_set(J, iter, 5, 0);
            gsl_matrix_set(J, iter, 6, 0);
        }
        return GSL_SUCCESS;
    }
    if (lambda<0)
    {
        for (int iter=0; iter<num_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 0);
            gsl_matrix_set(J, iter, 2, 0);
            gsl_matrix_set(J, iter, 3, 2e5*lambda);
            gsl_matrix_set(J, iter, 4, 0);
            gsl_matrix_set(J, iter, 5, 0);
            gsl_matrix_set(J, iter, 6, 0);
        }
        return GSL_SUCCESS;
    }
    
//    cout << q << endl;
//    for (int iterpara=0; iterpara<numOfPara; ++iterpara)
//    {
//        cout << gsl_vector_get(para, iterpara) << endl;
//    }
//    cout << endl;
    
    //Temperary variables used for acceleration.
    long double Dq2lambda=Dq2+lambda;
    long double vbsigma2=vbar/sigma/sigma;
    long double vb2sigma2=vbsigma2*vbar;
    long double logvbsigma2=log(vbsigma2);
    long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    long double vb2sigma3=vb2sigma2/sigma;
    
    long double paraISF[9]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2, vbar, cpsiz1, vb2sigma3};
    
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
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "dvisf" << endl;
    dvbarILT->NiLT_weeks(paraISF);
    //cout << "dzisf" << endl;
    dsigmaILT->NiLT_weeks(paraISF);
    //cout << "ddisf" << endl;
    dDILT->NiLT_weeks(paraISF);
    //cout << "dlisf" << endl;
    dlambdaILT->NiLT_weeks(paraISF);
    //cout << "ok" << endl;
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        long double t=tau[iter];
        double rtd=ILT->clenshaw(t);
        double dvbarrtd=dvbarILT->clenshaw(t);
        double dsigmartd=dsigmaILT->clenshaw(t);
        double dlambdartd=dlambdaILT->clenshaw(t);
        double dDrtd=dDILT->clenshaw(t);
        
        double expterm=exp(-Dq2*t);
        double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
        double yi=A*dA+B;
        
        gsl_matrix_set(J, iter, 0, A*(expterm-rtd)/yi );
        gsl_matrix_set(J, iter, 1, -A*alpha*dvbarrtd/yi );
        gsl_matrix_set(J, iter, 2, -A*alpha*dsigmartd/yi );
        gsl_matrix_set(J, iter, 3, -A*alpha*dlambdartd/yi );
        
        gsl_matrix_set(J, iter, 4, A*((1.0-alpha)*expterm*q*q*t-alpha*dDrtd)/yi );
        
        gsl_matrix_set(J, iter, 5, dA/yi );
        gsl_matrix_set(J, iter, 6, 1.0/yi );
        
        //Punishment terms
        if (alpha<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *alpha);
        }
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
        }
        if (D<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *D);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 5, gsl_matrix_get(J, iter, 5) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#else
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double* qArray=((dataStruct *)sdata)->q;
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    
    int* num_fitArray=((dataStruct *)sdata)->num_fit;
    
    long double q=qArray[0];
    int num_fit1=num_fitArray[0];
    int num_fit2=num_fit1+num_fitArray[1];
    
    //Get the parameters.
    long double alpha=gsl_vector_get(para, 0);
    long double vbar=gsl_vector_get(para, 1);
    long double sigma=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    long double A=gsl_vector_get(para, 5);
    long double B=gsl_vector_get(para, 6);
    
    long double Dq2=D*q*q;
    
    if (sigma<0)
    {
        for (int iter = 0; iter<num_fit2; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*sigma*sigma - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    if (vbar<0)
    {
        for (int iter = 0; iter<num_fit2; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*vbar*vbar - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    if (lambda<0)
    {
        for (int iter = 0; iter<num_fit2; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*lambda*lambda - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    
    //cout << q << endl;
    //for (int iterpara=0; iterpara<numOfPara; ++iterpara)
    //{
        //cout << gsl_vector_get(para, iterpara) << endl;
    //}
    //cout << endl;
    
    ////cout << lambda << endl;
    
    //Temperary variables used for acceleration.
    long double Dq2lambda=Dq2+lambda;
    long double vbsigma2=vbar/sigma/sigma;
    long double vb2sigma2=vbsigma2*vbar;
    long double logfactor=vb2sigma2*log(vbsigma2)-gsl_sf_lngamma(vb2sigma2);
    
    long double paraISF[6]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2};
    
    //Initialization of numerical inverse Laplace transformation solver
    int tid=omp_get_thread_num();
    ILT->cfun[tid].fun=ISFs;
    
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    //Set sigma and b in iLT solver
    csigma=-Dq2+0.5l;
    cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
    
    //Temperary variables used for acceleration.
    cb2=cb*2;
    csigmab=csigma-cb;
    
    //Calculate the coefficients of Laguerre polynomial series expansion.
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "ok" << endl;
    
    for (int iter = 0; iter<num_fit1; ++iter)
    {
        long double t=tau[iter];
        //Evaluate ISF at time t, the coefficients has been calculated.
        double rtd=ILT->clenshaw(t);
        double yi=log(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
        
        double result = yi - dataAry[iter];
        //Punishment terms, to make constrains in parameter space.
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
    
    q=qArray[1];
    
    //Get the parameters.
    A=gsl_vector_get(para, 7);
    B=gsl_vector_get(para, 8);
    
    Dq2=D*q*q;
    
    //cout << q << endl;
    //for (int iterpara=0; iterpara<numOfPara; ++iterpara)
    //{
        //cout << gsl_vector_get(para, iterpara) << endl;
    //}
    //cout << endl;
    
    ////cout << lambda << endl;
    
    //Temperary variables used for acceleration.
    Dq2lambda=Dq2+lambda;
    
    paraISF[1]=q;
    paraISF[2]=Dq2lambda;
    
    //Set sigma and b in iLT solver
    csigma=-Dq2+0.5l;
    cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
    
    //Temperary variables used for acceleration.
    cb2=cb*2;
    csigmab=csigma-cb;
    
    //Calculate the coefficients of Laguerre polynomial series expansion.
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "ok" << endl;
    
    for (int iter = num_fit1; iter<num_fit2; ++iter)
    {
        long double t=tau[iter];
        //Evaluate ISF at time t, the coefficients has been calculated.
        double rtd=ILT->clenshaw(t);
        double yi=log(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
        
        double result = yi - dataAry[iter];
        //Punishment terms, to make constrains in parameter space.
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
    double* qArray=((dataStruct *)sdata)->q;
    
    int* num_fitArray=((dataStruct *)sdata)->num_fit;
    
    long double q=qArray[0];
    int num_fit1=num_fitArray[0];
    int num_fit2=num_fit1+num_fitArray[1];
    
    long double alpha=gsl_vector_get(para, 0);
    long double vbar=gsl_vector_get(para, 1);
    long double sigma=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    long double A=gsl_vector_get(para, 5);
    long double B=gsl_vector_get(para, 6);
    
    long double Dq2=D*q*q;
    
    if (sigma<0)
    {
        for (int iter=0; iter<num_fit2; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 0);
            gsl_matrix_set(J, iter, 2, 2e5 *sigma);
            gsl_matrix_set(J, iter, 3, 0);
            gsl_matrix_set(J, iter, 4, 0);
            gsl_matrix_set(J, iter, 5, 0);
            gsl_matrix_set(J, iter, 6, 0);
        }
        return GSL_SUCCESS;
    }
    if (vbar<0)
    {
        for (int iter=0; iter<num_fit2; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 2e5 *vbar);
            gsl_matrix_set(J, iter, 2, 0);
            gsl_matrix_set(J, iter, 3, 0);
            gsl_matrix_set(J, iter, 4, 0);
            gsl_matrix_set(J, iter, 5, 0);
            gsl_matrix_set(J, iter, 6, 0);
        }
        return GSL_SUCCESS;
    }
    if (lambda<0)
    {
        for (int iter=0; iter<num_fit2; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 0);
            gsl_matrix_set(J, iter, 2, 0);
            gsl_matrix_set(J, iter, 3, 2e5*lambda);
            gsl_matrix_set(J, iter, 4, 0);
            gsl_matrix_set(J, iter, 5, 0);
            gsl_matrix_set(J, iter, 6, 0);
        }
        return GSL_SUCCESS;
    }
    
    //cout << q << endl;
    //for (int iterpara=0; iterpara<numOfPara; ++iterpara)
    //{
        //cout << gsl_vector_get(para, iterpara) << endl;
    //}
    //cout << endl;
    
    //Temperary variables used for acceleration.
    long double Dq2lambda=Dq2+lambda;
    long double vbsigma2=vbar/sigma/sigma;
    long double vb2sigma2=vbsigma2*vbar;
    long double logvbsigma2=log(vbsigma2);
    long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    long double vb2sigma3=vb2sigma2/sigma;
    
    long double paraISF[9]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2, vbar, cpsiz1, vb2sigma3};
    
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
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "dvisf" << endl;
    dvbarILT->NiLT_weeks(paraISF);
    //cout << "dzisf" << endl;
    dsigmaILT->NiLT_weeks(paraISF);
    //cout << "ddisf" << endl;
    dDILT->NiLT_weeks(paraISF);
    //cout << "dlisf" << endl;
    dlambdaILT->NiLT_weeks(paraISF);
    //cout << "ok" << endl;
    
    for (int iter=0; iter<num_fit1; ++iter)
    {
        //Temperary variables used for acceleration.
        long double t=tau[iter];
        double rtd=ILT->clenshaw(t);
        double dvbarrtd=dvbarILT->clenshaw(t);
        double dsigmartd=dsigmaILT->clenshaw(t);
        double dlambdartd=dlambdaILT->clenshaw(t);
        double dDrtd=dDILT->clenshaw(t);
        
        double expterm=exp(-Dq2*t);
        double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
        double yi=A*dA+B;
        
        gsl_matrix_set(J, iter, 0, A*(expterm-rtd)/yi );
        gsl_matrix_set(J, iter, 1, -A*alpha*dvbarrtd/yi );
        gsl_matrix_set(J, iter, 2, -A*alpha*dsigmartd/yi );
        gsl_matrix_set(J, iter, 3, -A*alpha*dlambdartd/yi );
        
        gsl_matrix_set(J, iter, 4, A*((1.0-alpha)*expterm*q*q*t-alpha*dDrtd)/yi );
        
        gsl_matrix_set(J, iter, 5, dA/yi );
        gsl_matrix_set(J, iter, 6, 1.0/yi );
        
        //Punishment terms
        if (alpha<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *alpha);
        }
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
        }
        if (D<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *D);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 5, gsl_matrix_get(J, iter, 5) + 2e5 *A);
        }
    }
    
    q=qArray[1];
    A=gsl_vector_get(para, 7);
    B=gsl_vector_get(para, 8);
    
    Dq2=D*q*q;
    Dq2lambda=Dq2+lambda;
    
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
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "dvisf" << endl;
    dvbarILT->NiLT_weeks(paraISF);
    //cout << "dzisf" << endl;
    dsigmaILT->NiLT_weeks(paraISF);
    //cout << "ddisf" << endl;
    dDILT->NiLT_weeks(paraISF);
    //cout << "dlisf" << endl;
    dlambdaILT->NiLT_weeks(paraISF);
    //cout << "ok" << endl;
    
    for (int iter=num_fit1; iter<num_fit2; ++iter)
    {
        //Temperary variables used for acceleration.
        long double t=tau[iter];
        double rtd=ILT->clenshaw(t);
        double dvbarrtd=dvbarILT->clenshaw(t);
        double dsigmartd=dsigmaILT->clenshaw(t);
        double dlambdartd=dlambdaILT->clenshaw(t);
        double dDrtd=dDILT->clenshaw(t);
        
        double expterm=exp(-Dq2*t);
        double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
        double yi=A*dA+B;
        
        gsl_matrix_set(J, iter, 0, A*(expterm-rtd)/yi );
        gsl_matrix_set(J, iter, 1, -A*alpha*dvbarrtd/yi );
        gsl_matrix_set(J, iter, 2, -A*alpha*dsigmartd/yi );
        gsl_matrix_set(J, iter, 3, -A*alpha*dlambdartd/yi );
        
        gsl_matrix_set(J, iter, 4, A*((1.0-alpha)*expterm*q*q*t-alpha*dDrtd)/yi );
        
        gsl_matrix_set(J, iter, 5, dA/yi );
        gsl_matrix_set(J, iter, 6, 1.0/yi );
        
        //Punishment terms
        if (alpha<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *alpha);
        }
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
        }
        if (D<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *D);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 5, gsl_matrix_get(J, iter, 5) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif

#endif
