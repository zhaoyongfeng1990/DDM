//
//  ISFSwimmer.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//
#include "ddm.h"

#ifdef ISFSW
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    const int num_fit=((dataStruct *)sdata)->num_fit;
    const int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    const double* tau=((dataStruct *)sdata)->tau;
    const double* qArray=((dataStruct *)sdata)->q;
    const double* dataAry=((dataStruct *)sdata)->data;
    
    //Get the parameters.
    const double alpha=gsl_vector_get(para, 0);
    const double D=gsl_vector_get(para, 1);
    const double v=gsl_vector_get(para, 2);
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        const double A=gsl_vector_get(para, 3+2*iterqc);
        const double B=gsl_vector_get(para, 4+2*iterqc);
        const double q=qArray[iterqc];
        
        for (int iter = 0; iter<num_fit; ++iter)
        {
            const int cidx=iter+iterqc*num_fit;
            const double t=tau[cidx];
            const double qvt=q*v*t;
            const double yi=(A*(1-exp(-D*q*q*t)*(1-alpha+alpha*sin(qvt)/qvt))+B);
            
            const double weight=1.0/sqrt(dataAry[cidx]);
            double result = (yi - dataAry[cidx])*weight;
            
            //Punishment terms, to make constrains in parameter space.
            if (alpha < 0)
            {
                result += 1e7 * alpha *alpha;
            }
            
            if (alpha>1)
            {
                result += 1e7*(alpha - 1)*(alpha - 1);
            }
            
            if (D<0)
            {
                result += 1e7*D*D;
            }
            
            if (v<0)
            {
                result += 1e7*v*v;
            }
            if (A<0)
            {
                result += 1e7*A*A;
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
    const int num_fit=((dataStruct *)sdata)->num_fit;
    const int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    const double* tau=((dataStruct *)sdata)->tau;
    const double* qArray=((dataStruct *)sdata)->q;
    const double* dataAry=((dataStruct *)sdata)->data;
    
    //Get the parameters
    const double alpha = gsl_vector_get(para, 0);
    const double D = gsl_vector_get(para, 1);
    const double v = gsl_vector_get(para, 2);
    
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        const double A=gsl_vector_get(para, 3+2*iterqc);
        const double B=gsl_vector_get(para, 4+2*iterqc);
        const double q=qArray[iterqc];
        for (int iter=0; iter<num_fit; ++iter)
        {
            const int cidx=iter+iterqc*num_fit;
            const double t=tau[cidx];
            //Temperary variables used for acceleration.
            const double qvt=q*v*t;
            const double sinqvt=sin(qvt);
            const double difexp=exp(-D*q*q*t);
            const double dyda=1-difexp*(1-alpha*alpha*sinqvt/qvt);
            const double weight=1.0/sqrt(dataAry[cidx]);
            //(log y)'=y'/y
            
            gsl_matrix_set(J, iter, 0, (A*difexp*(1-sinqvt/qvt))*weight );
            gsl_matrix_set(J, iter, 1, (A*q*q*t*(1-alpha+alpha*sinqvt/qvt)*difexp)*weight );
            
            gsl_matrix_set(J, iter, 2, A*difexp*alpha*(sinqvt/qvt-cos(qvt))/v*weight );
            
            gsl_matrix_set(J, iter, 3+2*iterqc, dyda*weight );
            gsl_matrix_set(J, iter, 4+2*iterqc, 1.0*weight);
            
            //Punishment terms, to make constrains in parameter space.
            if (alpha < 0)
            {
                gsl_matrix_set(J,cidx,0, gsl_matrix_get(J, cidx, 0)+ 2e7 *alpha);
            }
            
            if (alpha>1)
            {
                gsl_matrix_set(J, cidx, 0, gsl_matrix_get(J, cidx, 0) + 2e7 *(alpha-1));
            }
            
            if (D<0)
            {
                gsl_matrix_set(J, cidx, 1, gsl_matrix_get(J, cidx, 1) + 2e7 *D);
            }
            
            if (v<0)
            {
                gsl_matrix_set(J, cidx, 2, gsl_matrix_get(J, cidx, 2) + 2e7 *v);
            }
            
            if (A<0)
            {
                gsl_matrix_set(J, cidx, 3+2*iterqc, gsl_matrix_get(J, cidx, 3+2*iterqc) + 2e7 *A);
            }
        }
    }
    return GSL_SUCCESS;
}
#endif