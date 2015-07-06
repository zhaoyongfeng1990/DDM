//
//  ISFRunAndTumble3D.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"

#ifdef ISFRUNANDTUMBLE_3D

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    int num_fit=((dataStruct *)sdata)->num_fit;
    
    double v0=gsl_vector_get(para, 0);
    double lambda=gsl_vector_get(para, 1);
    double A=gsl_vector_get(para, 2);
    double B=gsl_vector_get(para, 3);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        double kv0=q*v0;
        double atanterm=atan(kv0/(tau[iter]+lambda));
        double yi=log(A*(1.0/tau[iter]-atanterm/(kv0-lambda*atanterm))+B/tau[iter]);
        
        double result = yi - dataAry[iter];
        //Punishment terms
        if (v0<0)
        {
            result += 1e5*v0*v0;
        }
        if (lambda<0)
        {
            result += 1e5*lambda*lambda;
        }
        if (A<0)
        {
            result += 1e5*A*A;
        }
        
        gsl_vector_set(y, iter, result);
    }
    
    return GSL_SUCCESS;
}

int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    int num_fit=((dataStruct *)sdata)->num_fit;
    
    double v0=gsl_vector_get(para, 0);
    double lambda=gsl_vector_get(para, 1);
    double A=gsl_vector_get(para, 2);
    double B=gsl_vector_get(para, 3);
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double kv0=q*v0;
        double lt=lambda+tau[iter];
        double atanterm=atan(kv0/lt);
        double lt2k2v02=lt*lt+kv0*kv0;
        double denominator2=kv0-lambda*atanterm;
        double denominator=lt2k2v02*denominator2*denominator2;
        double yi=A*(1.0/tau[iter]-atanterm/denominator2)+B/tau[iter];
        
        gsl_matrix_set(J, iter, 0, A*q*(lt2k2v02*atanterm-lt*kv0)/denominator/yi );
        
        gsl_matrix_set(J, iter, 1, A*(kv0*kv0-lt2k2v02*atanterm*atanterm)/denominator/yi );
        
        gsl_matrix_set(J, iter, 2, (1.0/tau[iter]-atanterm/denominator2)/yi );
        gsl_matrix_set(J, iter, 3, 1.0/tau[iter]/yi);
        
        //Punishment terms
        if (v0<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *v0);
        }
        if (lambda<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *lambda);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif
