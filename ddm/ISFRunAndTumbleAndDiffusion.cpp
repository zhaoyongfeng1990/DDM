//
//  ISFRunAndTumbleAndDiffusion.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "functions.h"

#ifdef ISFRunAndTumbleAndDiffusion
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters.
    double alpha=gsl_vector_get(para, 0);
    double v0=gsl_vector_get(para, 1);
    double lambda=gsl_vector_get(para, 2);
    double D=gsl_vector_get(para, 3);
    double A=gsl_vector_get(para, 4);
    double B=gsl_vector_get(para, 5);
    
    double kv0=q*v0;
    double dq2=D*q*q;
    for (int iter = 0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double atanterm=atan(kv0/(tau[iter]+lambda+dq2));
        double yi=log(A*(1.0/tau[iter]-(1.0-alpha)/(dq2+tau[iter])-alpha*atanterm/(kv0-lambda*atanterm))+B/tau[iter]);
        
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
    double q=((dataStruct *)sdata)->q;
    
    double alpha=gsl_vector_get(para, 0);
    double v0=gsl_vector_get(para, 1);
    double lambda=gsl_vector_get(para, 2);
    double D=gsl_vector_get(para, 3);
    double A=gsl_vector_get(para, 4);
    double B=gsl_vector_get(para, 5);
    
    double kv0=q*v0;
    double dq2=D*q*q;
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double tdq2=tau[iter]+dq2;
        double lt=lambda+tdq2;
        double atanterm=atan(kv0/lt);
        double lt2k2v02=lt*lt+kv0*kv0;
        double denominator2=kv0-lambda*atanterm;
        double denominator=lt2k2v02*denominator2*denominator2;
        double yi=A*(1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*atanterm/denominator2)+B/tau[iter];
        
        gsl_matrix_set(J, iter, 0, A*(1.0/tdq2-atanterm/denominator2)/yi );
        gsl_matrix_set(J, iter, 1, A*q*alpha*(lt2k2v02*atanterm-lt*kv0)/denominator/yi );
        
        gsl_matrix_set(J, iter, 2, A*alpha*(kv0*kv0-lt2k2v02*atanterm*atanterm)/denominator/yi );
        
        gsl_matrix_set(J, iter, 3, A*q*q*((1.0-alpha)/tdq2/tdq2+alpha*kv0*kv0/denominator)/yi );
        
        gsl_matrix_set(J, iter, 4, (1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*atanterm/denominator2)/yi );
        gsl_matrix_set(J, iter, 5, 1.0/tau[iter]/yi);
        
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
