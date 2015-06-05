//
//  ISFRunAndTumble2D.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "functions.h"

#ifdef ISFRUNANDTUMBLE
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    double lambda=gsl_vector_get(para, 0);
    double v = gsl_vector_get(para, 1);
    double A = gsl_vector_get(para, 2);
    double B = gsl_vector_get(para, 3);
    
    double kv0=q*v;
    for (int iter = 0; iter<num_fit; ++iter)
    {
        int n=2;
        double incrementeven = 1;
        double result = incrementeven*gsl_sf_bessel_J0(kv0*tau[iter]);
        double bracket = lambda*lambda*tau[iter] / kv0;
        double incrementodd = sqrt(bracket);
        result += incrementodd*gsl_sf_bessel_j0(kv0*tau[iter])*sqrt(kv0*tau[iter]);
        double err=0;
        //Calculation of ISF by using the infinite series expansion, small trick is used for numerical stability.
        do
        {
            int n2 = n / 2;
            incrementeven *= bracket / (n - 1);
            result += incrementeven*gsl_sf_bessel_Jn(n2, kv0*tau[iter]);
            
            incrementodd *= bracket / n;
            err = incrementodd*gsl_sf_bessel_jl(n2, kv0*tau[iter])*sqrt(kv0*tau[iter]);
            result += err;
            n+=2;
        } while (abs(err/result)>precision);
        result *= exp(-lambda*tau[iter]);
        result = A*(1 - result) + B;
        
        gsl_vector_set(y, iter, result - dataAry[iter]);
    }
    return GSL_SUCCESS;
}

int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau = ((dataStruct *)sdata)->tau;
    double q = ((dataStruct *)sdata)->q;
    
    double lambda = gsl_vector_get(para, 0);
    double v = gsl_vector_get(para, 1);
    double A = gsl_vector_get(para, 2);
    
    double kv0 = q*v;
    for (int iter = 0; iter<num_fit; ++iter)
    {
        double explt = exp(-lambda*tau[iter]);
        double bracket = lambda*lambda*tau[iter] / kv0;
        double relerr = 0;
        
        double incrementeven = 1;
        double j = gsl_sf_bessel_J0(kv0*tau[iter]);
        double f = j;
        double dfdl = tau[iter] * j;
        double dfdv = bracket*q*tau[iter]*j / 2.0; //j0 term
        
        //Calculation of Jacobian of ISF by using the infinite series expansion, small trick is used for numerical stability.
        
        double incrementodd = sqrt(bracket);
        j = incrementodd*gsl_sf_bessel_j0(kv0*tau[iter])*sqrt(kv0*tau[iter]);
        f += j;
        dfdl += (tau[iter] - 1.0 / lambda)*j;
        dfdv += j*q*tau[iter] / 4 * bracket;  //j1/2 term
        
        incrementeven *= bracket;
        j = gsl_sf_bessel_J1(kv0*tau[iter]);
        dfdv += q*tau[iter] * (bracket*bracket / 6.0 - 1)*j; //j1 term
        j *= incrementeven;
        f += j;
        dfdl += (tau[iter] - 2.0 / lambda)*j;
        
        incrementodd *= bracket / 2;
        j = incrementodd*gsl_sf_bessel_j1(kv0*tau[iter])*sqrt(kv0*tau[iter]);
        f += j;
        dfdl += j* (tau[iter] - 3.0 / lambda);
        dfdv += j*(lambda*lambda*tau[iter] * tau[iter] / 8.0 / v - q*q*v / lambda / lambda - 3.0 / 2.0 / v);
        
        dfdv -= lambda*q*tau[iter] * tau[iter] * gsl_sf_bessel_y0(kv0*tau[iter]); //y0 term
        int n = 4;
        do
        {
            int n2 = n / 2;
            incrementeven *= bracket / (n - 1);
            j = incrementeven*gsl_sf_bessel_Jn(n2, kv0*tau[iter]);
            f += j;
            dfdl += (tau[iter] - n / lambda)*j;
            dfdv += (lambda*lambda*tau[iter] * tau[iter] / 2.0/(n+1) / v - q*q*v*(n-1)/2.0 / lambda / lambda - n / 2.0 / v) *j;
            
            incrementodd *= bracket / n;
            j = incrementodd*gsl_sf_bessel_jl(n2, kv0*tau[iter])*sqrt(kv0*tau[iter]);
            f += j;
            relerr = abs(j / f);
            double ddfdl=j * (tau[iter] - (n + 1) / lambda);
            dfdl += ddfdl;
            if (ddfdl/dfdl>relerr)
            {
                relerr = ddfdl / dfdl;
            }
            double ddfdv = j* (lambda*lambda*tau[iter] * tau[iter] / 2.0 / (n + 2) / v - q*q*v*n / 2.0 / lambda / lambda - (n + 1) / 2.0 / v);
            dfdv += ddfdv;
            if (ddfdv / dfdv>relerr)
            {
                relerr = ddfdv / dfdv;
            }
            n += 2;
        } while (relerr>precision);
        f *= explt;
        dfdl *= explt;
        dfdv *= explt;
        
        gsl_matrix_set(J, iter, 0, A*dfdl);
        gsl_matrix_set(J, iter, 1, -A*dfdv);
        
        gsl_matrix_set(J, iter, 2, 1-f);
        gsl_matrix_set(J, iter, 3, 1);
    }
    return GSL_SUCCESS;
}

#endif