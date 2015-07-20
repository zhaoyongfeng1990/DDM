//
//  fitting.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_multifit_nlin.h>
#include <omp.h>

int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J);

int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole);
//Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.

int covar_rel_test(const gsl_matrix* J, const gsl_vector* x, double tol);
//Test fitting result using error estimation from covariance matrix, not reliable. tol is the relative tolerance of error.

void ddm::fitting()
{
    int cqsize=qsize-qIncreList[num_qCurve-1];
    int cnum_qCurve=num_qCurve;
    int cnum_fit=num_fit;
    int ctnum_fit=cnum_fit*num_qCurve;
    int cnumOfPara=numOfPara+2*num_qCurve;
    
    fittedPara=gsl_matrix_alloc(cqsize, cnumOfPara);
    //To store the fitting result and error.
    fitErr=gsl_matrix_alloc(cqsize, cnumOfPara);
    status = new int[cqsize];		//Record the status of fitting.
    
    const gsl_multifit_fdfsolver_type *solverType = gsl_multifit_fdfsolver_lmsder;
    //Using Levenberg-Marquardt algorithm as implemented in the scaled lmder routine in minpack. Jacobian is given.
    
    int progress=0;		//Indicator of progress.
    
#ifdef ISFRTD
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS), NILT3(OMP_NUM_THREADS), NILT4(OMP_NUM_THREADS);
#endif
    
#ifdef ISFRTDP
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS), NILT3(OMP_NUM_THREADS), NILT4(OMP_NUM_THREADS), NILT5(OMP_NUM_THREADS);
#endif
    
#pragma omp parallel for
    for (int iterq=0; iterq<cqsize; ++iterq)
    {
        //Get the selected data for fitting
        double* datafit=new double[ctnum_fit];
        double* qList=new double[cnum_qCurve];
        double* time=new double[ctnum_fit];
        for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
        {
            for (int iterf = 0; iterf < cnum_fit; ++iterf)
            {
                datafit[iterf+iterqc*cnum_fit]=log(gsl_matrix_get(datag, iterq+qIncreList[iterqc], iterf));		//Fitting in log scale.
                time[iterf+iterqc*cnum_fit]=tau[iterf];
            }
            qList[iterqc]=qabs[iterq]+qabs[iterq+qIncreList[iterqc]];
        }
        
        gsl_multifit_function_fdf fitfun;		//Function point.
        dataStruct sdata;		//GSL data structure
        
        sdata.data=datafit;
        sdata.tau=time;
        sdata.q=qList;
        sdata.num_fit=ctnum_fit;
        sdata.num_qCurve=cnum_qCurve;
        
#ifdef ISFRTD
        sdata.ISFILT=&NILT1;
        sdata.dvISFILT=&NILT2;
        sdata.dDISFILT=&NILT3;
        sdata.dlambdaISFILT=&NILT4;
#endif
        
#ifdef ISFRTDP
        sdata.ISFILT=&NILT1;
        sdata.dvbarISFILT=&NILT2;
        sdata.dsigmaISFILT=&NILT3;
        sdata.dDISFILT=&NILT4;
        sdata.dlambdaISFILT=&NILT5;
#endif
        
        //API
        fitfun.f=&ISFfun;
#ifdef NoJacobian
        fitfun.df=0;
        fitfun.fdf=0;
#else
        fitfun.df=&dISFfun;
        fitfun.fdf=&fdISFfun;
#endif
        fitfun.n=cnum_fit;
        fitfun.p=cnumOfPara;
        fitfun.params=&sdata;
        
        double localinipara[cnumOfPara];
        for (int iterp=0; iterp<numOfPara; ++iterp)
        {
            localinipara[iterp]=inipara[iterp];
        }
        //Estimation of A(q) and B(q)
        for (int iterqc=0; iterqc<num_qCurve; ++iterqc)
        {
            localinipara[numOfPara-1+2*iterqc] = gsl_matrix_get(datag, iterq+qIncreList[iterqc], 0);
            localinipara[numOfPara-2+2*iterqc] = gsl_matrix_get(datag, iterq+qIncreList[iterqc], cnum_fit-1)-localinipara[numOfPara-1+2*iterqc];
        }
        
        //Initiallization of the solver
        gsl_vector_view para=gsl_vector_view_array(localinipara, cnumOfPara);
        gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(solverType, cnum_fit, cnumOfPara);
        gsl_multifit_fdfsolver_set(solver, &fitfun, &para.vector);
        int iter=0;
        //gsl_vector* g=gsl_vector_alloc(numOfPara);
        //        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        //        {
        //            cout << gsl_vector_get(solver->x, iterpara) << '\n';
        //        }
        //        cout << '\n';
        
        do
        {
            gsl_multifit_fdfsolver_iterate(solver);		//Iterate one step.
            status[iterq] = norm0_rel_test(solver->dx, solver->x, 1e-10, 1e-10);  //Test the exiting condition
            
            //            for (int iterpara=0; iterpara<numOfPara; ++iterpara)
            //            {
            //                cout << gsl_vector_get(solver->x, iterpara) << '\n';
            //            }
            //            cout << '\n';
            //gsl_multifit_gradient(solver->J,solver->f, g);
            //status[iterq-1]=gsl_multifit_test_gradient(g, 1e-5);
            //			status[iterq - 1] = covar_rel_test(solver->J, solver->x, 1e-4);
            
            ++iter;
            if (iter>maxIter)
            {
                status[iterq]=GSL_EMAXITER;
            }
        } while (status[iterq] == GSL_CONTINUE);
        //gsl_vector_free(g);
        
        //Estimating the error.
        gsl_matrix* covar=gsl_matrix_alloc(cnumOfPara, cnumOfPara);
        gsl_multifit_covar(solver->J, 0.0, covar);
        for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)	//Record result.
        {
            gsl_matrix_set(fittedPara, iterq, iterpara, gsl_vector_get(solver->x, iterpara) );
            gsl_matrix_set(fitErr, iterq, iterpara, sqrt(gsl_matrix_get(covar, iterpara, iterpara)) );    //Not presice in log scale
        }
        gsl_matrix_free(covar);
        gsl_multifit_fdfsolver_free(solver);
        
        progress+=1;
        cout << "Fitted q=" << qabs[iterq] << " at iter=" << iter << ", " << 100.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << ", "<< gsl_strerror(status[iterq]) << "." << '\n';
        for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)
        {
            cout << gsl_matrix_get(fittedPara, iterq, iterpara) << '\n';
        }
        delete [] datafit;
        delete [] qList;
    }
}

#ifndef NoJacobian
int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J)
{
    ISFfun(para, sdata, y);
    dISFfun(para, sdata, J);
    
    return  GSL_SUCCESS;
}
#endif

//Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.
int ddm::norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole)
{
    int cnumOfPara=numOfPara+2*num_qCurve;
    //Check if some parameters become NAN or INF.
    for (int iter = 0; iter < cnumOfPara; ++iter)
    {
        if (!isfinite(x->data[iter]))
        {
            return GSL_EOVRFLW;
        }
    }
    for (int iter = 0; iter < cnumOfPara; ++iter)
    {
        //double relerr = abs(dx->data[iter]/* / x->data[iter]*/);
        if (abs(dx->data[iter]) > tol*abs(x->data[iter])+tole)
        {
            return GSL_CONTINUE;
        }
    }
    return GSL_SUCCESS;
}

//Test fitting result using error estimation from covariance matrix, not reliable. tol is the relative tolerance of error.
int covar_rel_test(const gsl_matrix* J, const gsl_vector* x, double tol)
{
    //Check if some parameters become NAN or INF.
    for (int iter = 0; iter < numOfPara; ++iter)
    {
        if (!isfinite(x->data[iter]))
        {
            return GSL_EOVRFLW;
        }
    }
    gsl_matrix* covar = gsl_matrix_alloc(numOfPara, numOfPara);
    gsl_multifit_covar(J, 0.0, covar);	//Get the covariance matrix.
    double fitErr[numOfPara];
    for (int iterpara = 0; iterpara<numOfPara; ++iterpara)	//err_i=\sqrt{c_{ii}}
    {
        fitErr[iterpara]=sqrt(gsl_matrix_get(covar, iterpara, iterpara));
    }
    gsl_matrix_free(covar);
    for (int iter = 0; iter < numOfPara; ++iter)
    {
        double relerr = abs(fitErr[iter] / x->data[iter]);
        if (relerr>tol)
        {
            return GSL_CONTINUE;
        }
    }
    return GSL_SUCCESS;
}