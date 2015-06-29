//
//  ISFRTDWithP.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"

#ifdef ISFRunAndTumbleAndDiffusionAndPv
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_integration.h>

double integrandFun(double t, void* params);
double dvbarFun(double t, void* params);
double dZFun(double t, void* params);
double dlambdaFun(double t, void* params);
double dDFun(double t, void* params);

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
double integrandFun(double t, void* params)
{
    double* para=(double*)params;
    double Z=para[1];
    double lambda=para[2];
    double q=para[3];
    
    double v=(1.-t)/t;
    
    double kv=q*v;
    double atanterm=atan(kv/para[4]);
    double Z1vbar=para[5];
    double vZ1vbar=v*Z1vbar;
    double constant=para[6];
    
    double exponent=Z*log(v)+constant-vZ1vbar;

    return atanterm/(kv-lambda*atanterm)*exp(exponent)/t/t;
}

double dvbarFun(double t, void* params)
{
    double* para=(double*)params;
    double vbar=para[0];
    double Z=para[1];
    double lambda=para[2];
    double q=para[3];
    
    double v=(1.-t)/t;
    
    double kv=q*v;
    double atanterm=atan(kv/para[4]);
    double Z1vbar=para[5];
    double vZ1vbar=v*Z1vbar;
    double constant=para[6];
    double factor=para[7];      //Z1/vbar/vbar
    
    double exponent=Z*log(v)+constant-vZ1vbar;

    return atanterm/(kv-lambda*atanterm)*exp(exponent)*(v-vbar)*factor/t/t;
}

double dZFun(double t, void* params)
{
    double* para=(double*)params;
    double vbar=para[0];
    double Z=para[1];
    double lambda=para[2];
    double q=para[3];
    
    double v=(1.-t)/t;
    
    double kv=q*v;
    double atanterm=atan(kv/para[4]);
    double Z1vbar=para[5];
    double vZ1vbar=v*Z1vbar;
    double psiz1=para[8];
    double constant=para[6];
    
    double exponent=Z*log(v)+constant-vZ1vbar;
    
    return atanterm/(kv-lambda*atanterm)*exp(exponent)*(1-v/vbar+log(vZ1vbar)-psiz1)/t/t;
}

double dlambdaFun(double t, void* params)
{
    double* para=(double*)params;
    double Z=para[1];
    double lambda=para[2];
    double q=para[3];
    
    double v=(1.-t)/t;
    
    double kv=q*v;
    double kv2=kv*kv;
    double ldq2s=para[4];
    double atanterm=atan(kv/ldq2s);
    double kvlatan=kv-lambda*atanterm;
    double Z1vbar=para[5];
    double vZ1vbar=v*Z1vbar;
    double ldq2s2=para[10];
    double kv2ldq2s2=kv2+ldq2s2;
    double constant=para[6];
    
    double exponent=Z*log(v)+constant-vZ1vbar;
    
    return exp(exponent)*(atanterm*atanterm*kv2ldq2s2-kv2)/kvlatan/kvlatan/t/t/kv2ldq2s2;
}

double dDFun(double t, void* params)
{
    double* para=(double*)params;
    double Z=para[1];
    double lambda=para[2];
    double q=para[3];
    
    double v=(1.-t)/t;
    
    double kv=q*v;
    double ldq2s=para[4];
    double atanterm=atan(kv/ldq2s);
    double kvlatan=kv-lambda*atanterm;
    double Z1vbar=para[5];
    double vZ1vbar=v*Z1vbar;
    double ldq2s2=para[10];
    double factor=para[9];     //pow(q,4)
    double constant=para[6];
    
    double exponent=Z*log(v)+constant-vZ1vbar;

    return -exp(exponent)*factor*v*v/(kv*kv+ldq2s2)/kvlatan/kvlatan/t/t;
}

int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters.
    double IntPara[11];
    double alpha=gsl_vector_get(para, 0);
    double vbar=gsl_vector_get(para, 1);
    double Z=gsl_vector_get(para, 2);
    double lambda=0; //gsl_vector_get(para, 3);
    double D=gsl_vector_get(para, 4);
    double A=gsl_vector_get(para, 5);
    double B=gsl_vector_get(para, 6);

    //cout << q << endl;
    //cout << alpha << endl;
    //cout << vbar << endl;
    //cout << Z << endl;
    //cout << lambda << endl;
    //cout << D << endl;
    //cout << A << endl;
    //cout << B << endl;
    
    //IntPara[0]=vbar;
    IntPara[1]=Z;
    IntPara[2]=lambda;
    
    double dq2=D*q*q;
    IntPara[3]=q;
    double Z1=Z+1;
    
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
    if (lambda<-dq2)
    {
        for (int iter = 0; iter<num_fit; ++iter)
        {
            gsl_vector_set(y, iter, 1e5*lambda*lambda - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    
    IntPara[5]=Z1/vbar;
    
    IntPara[6]=Z1*log(IntPara[5])-gsl_sf_lngamma(Z1);
    //IntPara[7]=Z1/vbar/vbar;
    
    //IntPara[8]=gsl_sf_psi(Z1);
    
    //IntPara[9]=pow(q,4);
    
    gsl_function pFun;
    pFun.function=&integrandFun;
    pFun.params=IntPara;

    gsl_integration_workspace* workspace=gsl_integration_workspace_alloc(10000);
    //gsl_integration_cquad_workspace* workspace=gsl_integration_cquad_workspace_alloc(1000000);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        IntPara[4]=tau[iter]+lambda+dq2;
        //IntPara[10]=IntPara[4]*IntPara[4];
        
        double pts[2]={0,1};
        double yi=0;
        double err=0;
        //size_t nevals;
        gsl_integration_qagp(&pFun, pts, 2, 0, 1e-5, 10000, workspace, &yi, &err);
        //gsl_integration_cquad(&pFun, 1e-5, 1-1e-5, 0, 1e-5, workspace, &yi, &err, &nevals);
        yi=log(A*(1.0/tau[iter]-(1.0-alpha)/(dq2+tau[iter])-alpha*yi)+B/tau[iter]);

        double result = yi - dataAry[iter];
        //Punishment terms, to make constrains in parameter space.
        if (vbar<0)
        {
            result += 1e6*vbar*vbar;
        }
        if (lambda<0)
        {
            result += 1e6*lambda*lambda;
        }
        if (D<0)
        {
            result += 1e6*D*D;
        }
        if (alpha<0)
        {
            result += 1e6*alpha*alpha;
        }
        if (alpha>1)
        {
            result += 1e6*(alpha-1)*(alpha-1);
        }
        if (A<0)
        {
            result += 1e6*A*A;
        }
        if (Z<0)
        {
            result += 1e6*Z*Z;
        }

        gsl_vector_set(y, iter, result);
    }

    gsl_integration_workspace_free(workspace);
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    double IntPara[11];
    double alpha=gsl_vector_get(para, 0);
    double vbar=gsl_vector_get(para, 1);
    double Z=gsl_vector_get(para, 2);
    double lambda=0; //gsl_vector_get(para, 3);
    double D=gsl_vector_get(para, 4);
    double A=gsl_vector_get(para, 5);
    double B=gsl_vector_get(para, 6);
    
    IntPara[0]=vbar;
    IntPara[1]=Z;
    IntPara[2]=lambda;
    
    double dq2=D*q*q;
    IntPara[3]=q;
    double Z1=Z+1;
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
    if (lambda<-dq2)
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
    
    IntPara[5]=Z1/vbar;
    
    IntPara[6]=Z1*log(IntPara[5])-gsl_sf_lngamma(Z1);
    IntPara[7]=Z1/vbar/vbar;

    IntPara[8]=gsl_sf_psi(Z1);

    IntPara[9]=pow(q,4);
    
    double pts[2]={0,1};
    
    gsl_function pFun;
    pFun.function=&integrandFun;
    pFun.params=IntPara;
    
    gsl_function pdvbarFun;
    pdvbarFun.function=&dvbarFun;
    pdvbarFun.params=IntPara;
    
    gsl_function pdZFun;
    pdZFun.function=&dZFun;
    pdZFun.params=IntPara;
    
    gsl_function pdlambdaFun;
    pdlambdaFun.function=&dlambdaFun;
    pdlambdaFun.params=IntPara;
    
    gsl_function pdDFun;
    pdDFun.function=&dDFun;
    pdDFun.params=IntPara;
    
    gsl_integration_workspace* workspace=gsl_integration_workspace_alloc(10000);
    //gsl_integration_cquad_workspace* workspace=gsl_integration_cquad_workspace_alloc(1000000);
    for (int iter=0; iter<num_fit; ++iter)
    {
        IntPara[4]=tau[iter]+lambda+dq2;
        IntPara[10]=IntPara[4]*IntPara[4];
        
        double inteResult=0;
        double err=0;
        //size_t nevals;
        gsl_integration_qagp(&pFun, pts, 2, 0, 1e-5, 10000, workspace, &inteResult, &err);
        //gsl_integration_cquad(&pFun, 1e-5, 1-1e-5, 0, 1e-5, workspace, &inteResult, &err, &nevals);

        double dvbar=0;
        gsl_integration_qagp(&pdvbarFun, pts, 2, 0, 1e-5, 10000, workspace, &dvbar, &err);
        //gsl_integration_cquad(&pdvbarFun, 1e-5, 1-1e-5, 0, 1e-5, workspace, &dvbar, &err, &nevals);

        double dZ=0;
        gsl_integration_qagp(&pdZFun, pts, 2, 0, 1e-5, 10000, workspace, &dZ, &err);
        //gsl_integration_cquad(&pdZFun, 1e-5, 1-1e-5, 0, 1e-5, workspace, &dZ, &err, &nevals);

        //double dlambda=0;
        //gsl_integration_qagp(&pdlambdaFun, pts, 2, 0, 1e-5, 10000, workspace, &dlambda, &err);
        //gsl_integration_cquad(&pdlambdaFun, 0, 1, 0, 1e-5, workspace, &dlambda, &err, &nevals);

        double dD=0;
        gsl_integration_qagp(&pdDFun, pts, 2, 0, 1e-5, 10000, workspace, &dD, &err);
        //gsl_integration_cquad(&pdDFun, 1e-5, 1-1e-5, 0, 1e-5, workspace, &dD, &err, &nevals);

        //Temperary variables used for acceleration.
        double tdq2=tau[iter]+dq2;
        double dydA=(1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*inteResult);
        double yi=A*dydA+B/tau[iter];
        
        gsl_matrix_set(J, iter, 0, A*(1.0/tdq2-inteResult)/yi );
        gsl_matrix_set(J, iter, 1, -A*alpha*dvbar/yi );
        gsl_matrix_set(J, iter, 2, -A*alpha*dZ/yi );
        
        gsl_matrix_set(J, iter, 3, 0/*-A*alpha*dlambda/yi */);
        
        gsl_matrix_set(J, iter, 4, A*((1.0-alpha)*q*q/tdq2/tdq2-alpha*dD)/yi );
        
        gsl_matrix_set(J, iter, 5, dydA/yi );
        gsl_matrix_set(J, iter, 6, 1.0/tau[iter]/yi);
        
        //Punishment terms
        if (alpha<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0)+ 2e6 *alpha);
        }
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0)+ 2e6 *(alpha-1));
        }
        if (vbar<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1)+ 2e6 *vbar);
        }
        if (Z<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2)+ 2e6 *Z);
        }
        if (lambda<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3)+ 2e6 *lambda);
        }
        if (D<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4)+ 2e6 *D);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 5, gsl_matrix_get(J, iter, 5)+ 2e6 *A);
        }
    }
    
    gsl_integration_workspace_free(workspace);
    return GSL_SUCCESS;
}
#endif
