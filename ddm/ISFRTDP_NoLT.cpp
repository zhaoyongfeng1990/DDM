//
//  ISFRTDP_NoLT.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/26.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>

#ifdef ISFRTDPNoLT

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
    long double Z=para[3];
    long double logfactor=para[4];
    long double Z1vbar=para[5];
    
    long double qv=q*v;
    long double pv=exp(Z*log(v)+logfactor-v*Z1vbar);
    
    cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return atanterm/(qv-lambda*atanterm)*pv;
}

cpx dvbarISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double Z=para[3];
    long double logfactor=para[4];
    long double Z1vbar=para[5];
    long double vb=para[6];
    
    long double qv=q*v;
    long double pv=exp(Z*log(v)+logfactor-v*Z1vbar);
    
    cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return atanterm/(qv-lambda*atanterm)*pv*(v-vb)*Z1vbar/vb;
}

cpx dZISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double Z=para[3];
    long double logfactor=para[4];
    long double Z1vbar=para[5];
    long double vb=para[6];
    long double psiz1=para[7];
    
    long double vZ1vbar=v*Z1vbar;
    long double qv=q*v;
    
    long double pv=exp(Z*log(v)+logfactor-vZ1vbar);
    
    cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return atanterm/(qv-lambda*atanterm)*pv*(1.0l-v/vb+log(vZ1vbar)-psiz1);
}

cpx dDISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    long double lambda=para[0];
    long double q=para[1];
    long double Dq2lambda=para[2];
    long double Z=para[3];
    long double logfactor=para[4];
    long double Z1vbar=para[5];
    
    long double qv=q*v;
    long double qv2=qv*qv;
    long double pv=exp(Z*log(v)+logfactor-v*Z1vbar);
    
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
    long double Z=para[3];
    long double logfactor=para[4];
    long double Z1vbar=para[5];
    
    long double qv=q*v;
    long double qv2=qv*qv;
    long double pv=exp(Z*log(v)+logfactor-v*Z1vbar);
    
    cpx Dq2lambdas=Dq2lambda+s;
    cpx atanterm=atan(qv/Dq2lambdas);
    return ((qv2+Dq2lambdas*Dq2lambdas)*atanterm*atanterm-qv2)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm)*pv;
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
    long double vbar=gsl_vector_get(para, 1);
    long double Z=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    long double A=gsl_vector_get(para, 5);
    long double B=gsl_vector_get(para, 6);
    
    long double Dq2=D*q*q;
    
    if (Z<0)
    {
        for (int iter = 0; iter<num_fit; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*Z*Z - dataAry[iter]);
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
    
    //cout << q << endl;
    for (int iterpara=0; iterpara<numOfPara; ++iterpara)
    {
        //cout << gsl_vector_get(para, iterpara) << endl;
    }
    //cout << endl;
    
    ////cout << lambda << endl;
    
    //Temperary variables used for acceleration.
    long double Dq2lambda=Dq2+lambda;
    long double Z1=Z+1.0l;
    long double Z1vbar=Z1/vbar;
    long double logfactor=Z1*log(Z1vbar)-gsl_sf_lngamma(Z1);
    
    long double paraISF[6]={lambda, q, Dq2lambda, Z, logfactor, Z1vbar};
    
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
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        long double t=tau[iter];
        //Evaluate ISF at time t, the coefficients has been calculated.
        double rtd=ILT->clenshaw(t);
        double yi=log(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
        
        double result = yi - dataAry[iter];
        //Punishment terms, to make constrains in parameter space.
        if (vbar<0)
        {
            result += 1e5*vbar*vbar;
        }
        if (Z<0)
        {
            result += 1e5*Z*Z;
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
    long double vbar=gsl_vector_get(para, 1);
    long double Z=gsl_vector_get(para, 2);
    long double lambda=gsl_vector_get(para, 3);
    long double D=gsl_vector_get(para, 4);
    long double A=gsl_vector_get(para, 5);
    long double B=gsl_vector_get(para, 6);
    
    long double Dq2=D*q*q;
    
    if (Z<0)
    {
        for (int iter=0; iter<num_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 0);
            gsl_matrix_set(J, iter, 1, 0);
            gsl_matrix_set(J, iter, 2, 2e5 *Z);
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
    
    //cout << q << endl;
    for (int iterpara=0; iterpara<numOfPara; ++iterpara)
    {
        //cout << gsl_vector_get(para, iterpara) << endl;
    }
    //cout << endl;
    
    //Temperary variables used for acceleration.
    long double Dq2lambda=Dq2+lambda;
    long double Z1=Z+1.0l;
    long double Z1vbar=Z1/vbar;
    long double logfactor=Z1*log(Z1vbar)-gsl_sf_lngamma(Z1);
    long double psiz1=gsl_sf_psi(Z1);
    
    long double paraISF[8]={lambda, q, Dq2lambda, Z, logfactor, Z1vbar, vbar, psiz1};
    
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    NILT* dvbarILT=((dataStruct *)sdata)->dvbarISFILT;
    NILT* dZILT=((dataStruct *)sdata)->dZISFILT;
    NILT* dlambdaILT=((dataStruct *)sdata)->dlambdaISFILT;
    NILT* dDILT=((dataStruct *)sdata)->dDISFILT;
    
    int tid=omp_get_thread_num();
    ILT->cfun[tid].fun=ISFs;
    dvbarILT->cfun[tid].fun=dvbarISFs;
    dZILT->cfun[tid].fun=dZISFs;
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
    dZILT->sigma[tid]=csigma;
    dlambdaILT->sigma[tid]=csigma;
    dDILT->sigma[tid]=csigma;
    
    dvbarILT->b[tid]=cb;
    dZILT->b[tid]=cb;
    dlambdaILT->b[tid]=cb;
    dDILT->b[tid]=cb;
    
    dvbarILT->b2[tid]=cb2;
    dZILT->b2[tid]=cb2;
    dlambdaILT->b2[tid]=cb2;
    dDILT->b2[tid]=cb2;
    
    dvbarILT->sigmab[tid]=csigmab;
    dZILT->sigmab[tid]=csigmab;
    dlambdaILT->sigmab[tid]=csigmab;
    dDILT->sigmab[tid]=csigmab;
    
    //cout << "isf" << endl;
    ILT->NiLT_weeks(paraISF);
    //cout << "dvisf" << endl;
    dvbarILT->NiLT_weeks(paraISF);
    //cout << "dzisf" << endl;
    dZILT->NiLT_weeks(paraISF);
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
        double dZrtd=dZILT->clenshaw(t);
        double dlambdartd=dlambdaILT->clenshaw(t);
        double dDrtd=dDILT->clenshaw(t);
        
        double expterm=exp(-Dq2*t);
        double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
        double yi=A*dA+B;
        
        gsl_matrix_set(J, iter, 0, A*(expterm-rtd)/yi );
        gsl_matrix_set(J, iter, 1, -A*alpha*dvbarrtd/yi );
        gsl_matrix_set(J, iter, 2, -A*alpha*dZrtd/yi );
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
        if (vbar<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *vbar);
        }
        if (Z<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *Z);
        }
        if (lambda<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *lambda);
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
