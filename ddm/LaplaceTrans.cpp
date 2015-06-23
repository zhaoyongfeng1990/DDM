//
//  LaplaceTrans.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <gsl/gsl_blas.h>
#include <omp.h>

void ddm::LaplaceTrans()
{
    gsl_matrix* transM=gsl_matrix_alloc(numOfDiff, num_fit);		//Transformation matrix
    
#pragma omp parallel for
    for (int itert = 0; itert < numOfDiff; ++itert)
    {
        for (int iters = 0; iters < num_fit; ++iters)
        {
            double tau=(itert+1)*dt;
            s[iters]=smin+iters*ds;             //s is sampled in linear scale
            //s[iters]=exp(1+0.0005*iters);		//s is sampled in log scale
            gsl_matrix_set(transM, itert, iters, dt*exp(-tau*s[iters]));		//\Delta t=0.01, numerical integration here.
        }
    }
    
    gsl_vector_view lastRow=gsl_matrix_row(transM, numOfDiff-1);
    gsl_vector_scale(&lastRow.vector, 0.5);
    
    ldatag=gsl_matrix_alloc(qsize, num_fit);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, datag, transM, 0, ldatag);		//LT is performed by matrix multiplication. Call BLAS to calculate matrix multiplication.
    gsl_matrix_add_constant(ldatag, 0.5*dt);
    gsl_matrix_free(transM);
}