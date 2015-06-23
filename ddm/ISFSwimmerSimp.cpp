//
//  ISFSwimmer.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"

#ifdef ISFSWIMMERSIMPLE
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
    double A=gsl_vector_get(para, 3);
    double B=gsl_vector_get(para, 4);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        double qvt=q*v*tau[iter];
        double yi=log(A*(1-exp(-D*q*q*tau[iter])*(1-alpha+alpha*sin(qvt)/qvt))+B);
        
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
    double A = gsl_vector_get(para, 3);
    double B = gsl_vector_get(para, 4);
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double qvt=q*v*tau[iter];
        double sinqvt=sin(qvt);
        double difexp=exp(-D*q*q*tau[iter]);
        double dyda=1-difexp*(1-alpha*alpha*sinqvt/qvt);
        double yi = A*dyda + B;
        //(log y)'=y'/y
        
        gsl_matrix_set(J, iter, 0, (A*difexp*(1-sinqvt/qvt))/yi );
        gsl_matrix_set(J, iter, 1, (A*q*q*tau[iter]*(1-alpha+alpha*sinqvt/qvt)*difexp)/yi );
        
        gsl_matrix_set(J, iter, 2, A*difexp*alpha*(sinqvt/qvt-cos(qvt))/v/yi );
        
        gsl_matrix_set(J, iter, 3, dyda/yi );
        gsl_matrix_set(J, iter, 4, 1/yi);
        
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
        
        if (A<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif