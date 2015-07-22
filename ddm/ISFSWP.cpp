//
//  ISFSwimmer.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"
#include <iostream>
#ifdef ISFSWP

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    int num_fit=((dataStruct *)sdata)->num_fit;
    int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    double* tau=((dataStruct *)sdata)->tau;
    double* qArray=((dataStruct *)sdata)->q;
    double* dataAry=((dataStruct *)sdata)->data;
    
    //Get the parameters.
    double alpha=gsl_vector_get(para, 0);
    double D=gsl_vector_get(para, 1);
    double v=gsl_vector_get(para, 2);
    double Z=gsl_vector_get(para, 3);
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        double A=gsl_vector_get(para, 4+2*iterqc);
        double B=gsl_vector_get(para, 5+2*iterqc);
        double q=qArray[iterqc];
        for (int iter = 0; iter<num_fit; ++iter)
        {
            int cidx=iter+iterqc*num_fit;
            double t=tau[cidx];
            double Gamma=q*t*v/(Z+1);
            double yi=log(A*(1-exp(-D*q*q*t)*(1-alpha+alpha/Z/Gamma*sin(Z*atan(Gamma))/pow(1+Gamma*Gamma, Z/2.0)))+B);
            
            double weight=sqrt(exp(dataAry[cidx])/500/500);
            double result = (yi - dataAry[cidx])*weight;
            
            //Punishment terms, to make constrains in parameter space.
            if (alpha < 0)
            {
                result += 1e5 * alpha *alpha;
            }
            
            if (alpha>1)
            {
                result += 1e5*(alpha - 1)*(alpha - 1);
            }
            
            if (D<0)
            {
                result += 1e5*D*D;
            }
            
            if (v<0)
            {
                result += 1e5*v*v;
            }
            
            if (Z<0)
            {
                result += 1e5*Z*Z;
            }
            if (A<0)
            {
                result += 1e5*A*A;
            }
            
            gsl_vector_set(y, cidx, result);
        }
    }
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    gsl_matrix_set_zero(J);
    int num_fit=((dataStruct *)sdata)->num_fit;
    int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    double* tau=((dataStruct *)sdata)->tau;
    double* qArray=((dataStruct *)sdata)->q;
    double* dataAry=((dataStruct *)sdata)->data;
    
    //Get the parameters
    double alpha = gsl_vector_get(para, 0);
    double D = gsl_vector_get(para, 1);
    double v = gsl_vector_get(para, 2);
    double Z = gsl_vector_get(para, 3);
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        double A=gsl_vector_get(para, 4+2*iterqc);
        double B=gsl_vector_get(para, 5+2*iterqc);
        double q=qArray[iterqc];
        for (int iter=0; iter<num_fit; ++iter)
        {
            int cidx=iter+iterqc*num_fit;
            double t=tau[cidx];
            //Temperary variables used for acceleration.
            double Gamma=q*t*v/(Z+1);
            double difexp=exp(-D*q*q*t);
            double zatang=Z*atan(Gamma);
            double sinzg=sin(zatang);
            double coszg=cos(zatang);
            double powg=pow(1+Gamma*Gamma, Z/2.0);
            double qqttvv=q*q*t*t*v*v;
            double zzqqttvv = (1 + Z)*(1 + Z) + qqttvv;
            double dydA=(1 - difexp*(1 - alpha + alpha / Z / Gamma*sinzg / powg));
            double yi =A*dydA + B;
            //(log y)'=y'/y
            double weight=sqrt(exp(dataAry[cidx])/500/500);
            
            //cout << t << ": " << weight << '\n';
            
            gsl_matrix_set(J, cidx, 0, (A*difexp*(1-sinzg/Z/Gamma/powg))/yi*weight );
            gsl_matrix_set(J, cidx, 1, (A*q*q*t*(1-alpha+alpha*sinzg/powg/Z/Gamma)*difexp)/yi*weight );
            
            gsl_matrix_set(J, cidx, 2, ((alpha*(1+Z)*A*difexp*((1+Z+qqttvv)*sinzg-v*Z*q*t*coszg) )/powg/(v*Z*Gamma*zzqqttvv) )/yi*weight );
            
            gsl_matrix_set(J, cidx, 3, ((A*alpha*difexp*( (v*Z*Z*q*t-zzqqttvv*zatang)*2*coszg+(2*(1+Z-(Z-1)*qqttvv)+Z*zzqqttvv*log(1+Gamma*Gamma))*sinzg))/powg/Gamma/2/Z/Z/zzqqttvv)/yi*weight );
            
            gsl_matrix_set(J, cidx, 4+2*iterqc, dydA/yi*weight );
            gsl_matrix_set(J, cidx, 5+2*iterqc, 1/yi*weight );
            
            //Punishment terms, to make constrains in parameter space.
            if (alpha < 0)
            {
                gsl_matrix_set(J,cidx,0, gsl_matrix_get(J, cidx, 0)+ 2e5 *alpha);
            }
            
            if (alpha>1)
            {
                gsl_matrix_set(J, cidx, 0, gsl_matrix_get(J, cidx, 0) + 2e5 *(alpha-1));
            }
            
            if (D<0)
            {
                gsl_matrix_set(J, cidx, 1, gsl_matrix_get(J, cidx, 1) + 2e5 *D);
            }
            
            if (v<0)
            {
                gsl_matrix_set(J, cidx, 2, gsl_matrix_get(J, cidx, 2) + 2e5 *v);
            }
            
            if (Z<0)
            {
                gsl_matrix_set(J, cidx, 3, gsl_matrix_get(J, cidx, 3) + 2e5 *Z);
            }
            
            if (A<0)
            {
                gsl_matrix_set(J, cidx, 4+2*iterqc, gsl_matrix_get(J, cidx, 4+2*iterqc) + 2e5 *A);
            }
        }
    }
    return GSL_SUCCESS;
}
#endif