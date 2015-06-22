//
//  ISFRunAndTumbleAndDiffusionNoLT.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/18.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "functions.h"

#ifdef ISFRunAndTumbleAndDiffusionNoLT
#include <vector>
using namespace std;

typedef complex<long double> cpx;

cpx ISFs(cpx s, long double* para)
{
    long double qv=para[0];
    long double lambda=para[1];
    long double Dq2=para[2];
    long double q=para[3];
    cpx atanterm=atan(qv/(Dq2+s+lambda));
    return atanterm/(qv-lambda*atanterm);
}

void NiLT_weeks(cpx (*fun)(cpx, long double*), long double* para)
{
    for (int iter=0; iter<M; ++iter)
    {
        cpx itheta={0.0l,pi*(2.0l*iter+1.0l)/M};
        cpx expterm=exp(itheta)-1.0l;
        cpx s=sigma-b*(expterm+2.0l)/expterm;
        cpx r=-b2*fun(s,para)/expterm;
        iLTStruct.fftwIn[iter][0]=real(r);
        iLTStruct.fftwIn[iter][1]=imag(r);
    }
    fftwl_execute(iLTStruct.integration);
    for (int iter=0; iter<M; ++iter)
    {
        cpx temp={iLTStruct.fftwOut[iter][0],iLTStruct.fftwOut[iter][1]};
        cpx itheta={0,-iter*pi/M};
        temp*=exp(itheta);
        iLTStruct.CoeA[iter]=real(temp)/M;
    }
}

double clenshaw(long double t)
{
    long double y2=0.0l;
    long double y1=iLTStruct.CoeA[M-1];
    long double y0=0.0l;
    for (int k=M-1; k>0; --k)
    {
        long double lk=(long double)k;
        y0=(2.0l*lk-1.0l-b2*t)/lk*y1-lk/(lk+1.0l)*y2+iLTStruct.CoeA[k-1];
        y2=y1;
        y1=y0;
    }
    return exp(sigmab*t)*y0;
}

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    long double q=((dataStruct *)sdata)->q;
    
    //Get the parameters.
    long double alpha=0.8;//gsl_vector_get(para, 0);
    long double v0=13;//gsl_vector_get(para, 1);
    long double lambda=gsl_vector_get(para, 2);
    long double D=0.4;//gsl_vector_get(para, 3);
    long double A=1; //gsl_vector_get(para, 4);
    long double B=0; //gsl_vector_get(para, 5);
    
    //cout << lambda << endl;

    long double kv0=q*v0;
    long double Dq2=D*q*q;
    long double paraISF[4]={kv0, lambda, Dq2, q};
    NiLT_weeks(ISFs, paraISF);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        long double t=tau[iter];
        double rtd=clenshaw(t);
        //Temperary variables used for acceleration.
        double yi=log(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
        
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
//int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
//{
//    double* tau=((dataStruct *)sdata)->tau;
//    double q=((dataStruct *)sdata)->q;
//    
//    double alpha=0.8;//gsl_vector_get(para, 0);
//    double v0=13;//gsl_vector_get(para, 1);
//    double lambda=gsl_vector_get(para, 2);
//    double D=0.4;//gsl_vector_get(para, 3);
//    double A=1; //gsl_vector_get(para, 4);
//    double B=0; //gsl_vector_get(para, 5);
//    
//    double kv0=q*v0;
//    double dq2=D*q*q;
//    for (int iter=0; iter<num_fit; ++iter)
//    {
//        //Temperary variables used for acceleration.
//        double tdq2=tau[iter]+dq2;
//        double lt=lambda+tdq2;
//        double atanterm=atan(kv0/lt);
//        double lt2k2v02=lt*lt+kv0*kv0;
//        double denominator2=kv0-lambda*atanterm;
//        double denominator=lt2k2v02*denominator2*denominator2;
//        double yi=A*(1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*atanterm/denominator2)+B/tau[iter];
//        
//        gsl_matrix_set(J, iter, 0, 0/*A*(1.0/tdq2-atanterm/denominator2)/yi*/ );
//        gsl_matrix_set(J, iter, 1, 0/*A*q*alpha*(lt2k2v02*atanterm-lt*kv0)/denominator/yi*/ );
//        
//        gsl_matrix_set(J, iter, 2, A*alpha*(kv0*kv0-lt2k2v02*atanterm*atanterm)/denominator/yi );
//        
//        gsl_matrix_set(J, iter, 3, 0/*A*q*q*((1.0-alpha)/tdq2/tdq2+alpha*kv0*kv0/denominator)/yi*/ );
//        
//        gsl_matrix_set(J, iter, 4, 0 /*(1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*atanterm/denominator2)/yi */);
//        gsl_matrix_set(J, iter, 5, 0 /*1.0/tau[iter]/yi*/ );
//        
//        //Punishment terms
//        if (alpha<0)
//        {
//            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *alpha);
//        }
//        if (alpha>1)
//        {
//            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
//        }
//        if (v0<0)
//        {
//            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *v0);
//        }
//        if (lambda<0)
//        {
//            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *lambda);
//        }
//        if (D<0)
//        {
//            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *D);
//        }
//        if (A<0)
//        {
//            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *A);
//        }
//    }
//    
//    return GSL_SUCCESS;
//}
#endif
