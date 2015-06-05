//
//  ISFSwimmer.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "functions.h"

#ifdef ISFSWIMMER
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters.
    double alpha=gsl_vector_get(para, 0);
    double D=gsl_vector_get(para, 1);
    double v=gsl_vector_get(para, 2);
    double Z=gsl_vector_get(para, 3);
    double A=gsl_vector_get(para, 4);
    double B=gsl_vector_get(para, 5);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        double Gamma=q*tau[iter]*v/(Z+1);
        double yi=log(A*(1-exp(-D*q*q*tau[iter])*(1-alpha+alpha/Z/Gamma*sin(Z*atan(Gamma))/pow(1+Gamma*Gamma, Z/2.0)))+B);
        
        double result = yi - dataAry[iter];
        
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
        
        gsl_vector_set(y, iter, result);
    }
    
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters
    double alpha = gsl_vector_get(para, 0);
    double D = gsl_vector_get(para, 1);
    double v = gsl_vector_get(para, 2);
    double Z = gsl_vector_get(para, 3);
    double A = gsl_vector_get(para, 4);
    double B = gsl_vector_get(para, 5);
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double Gamma=q*tau[iter]*v/(Z+1);
        double difexp=exp(-D*q*q*tau[iter]);
        double zatang=Z*atan(Gamma);
        double sinzg=sin(zatang);
        double coszg=cos(zatang);
        double powg=pow(1+Gamma*Gamma, Z/2.0);
        double qqttvv=q*q*tau[iter]*tau[iter]*v*v;
        double zzqqttvv = (1 + Z)*(1 + Z) + qqttvv;
        double yi = A*(1 - difexp*(1 - alpha + alpha / Z / Gamma*sinzg / powg)) + B;
        //(log y)'=y'/y
        
        gsl_matrix_set(J, iter, 0, (A*difexp*(1-sinzg/Z/Gamma/powg))/yi );
        gsl_matrix_set(J, iter, 1, (A*q*q*tau[iter]*(1-alpha+alpha*sinzg/powg/Z/Gamma)*difexp)/yi );
        
        gsl_matrix_set(J, iter, 2, ((alpha*(1+Z)*A*difexp*((1+Z+qqttvv)*sinzg-v*Z*q*tau[iter]*coszg) )/powg/(v*Z*Gamma*zzqqttvv) )/yi );
        
        gsl_matrix_set(J, iter, 3, ((A*alpha*difexp*( (v*Z*Z*q*tau[iter]-zzqqttvv*zatang)*2*coszg+(2*(1+Z-(Z-1)*qqttvv)+Z*zzqqttvv*log(1+Gamma*Gamma))*sinzg))/powg/Gamma/2/Z/Z/zzqqttvv)/yi );
        
        gsl_matrix_set(J, iter, 4, (1-difexp*(1-alpha+alpha/Z/Gamma*sinzg/powg) )/yi );
        gsl_matrix_set(J, iter, 5, 1/yi);
        
        //Punishment terms, to make constrains in parameter space.
        if (alpha < 0)
        {
            gsl_matrix_set(J,iter,0, gsl_matrix_get(J, iter, 0)+ 2e5 *alpha);
        }
        
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
        }
        
        if (D<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *D);
        }
        
        if (v<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *v);
        }
        
        if (Z<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *Z);
        }
        
        if (A<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif