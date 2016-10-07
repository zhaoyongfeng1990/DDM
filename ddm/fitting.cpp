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

#ifdef ISFRTDPfix
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#endif

#ifdef ISFRTDPTTfix
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#endif

#include <omp.h>

int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J);

//Test fitting result using error estimation from covariance matrix, not reliable. tol is the relative tolerance of error.
int covar_rel_test(const gsl_matrix* J, const gsl_vector* x, double tol);

//Fitting. Allow fitting multiple q curves simultaneously to decrease the chance of converging to local minimum.
void ddm::fitting()
{
    int cnum_fit=num_fit;
    int ctimeWindow=timeWindow;
    //Find the truncation time if time window is set
    for (int itert=0; itert<num_fit; ++itert)
    {
        if (tau[itert]>ctimeWindow)
        {
            cnum_fit=itert;
            break;
        }
    }
    
    //Local variables
    int cqsize=qsize-qIncreList[num_qCurve-1];  //number of fitting result
    int cnum_qCurve=num_qCurve;
    int ctnum_fit=cnum_fit*num_qCurve;
    int cnumOfPara=numOfPara+2*num_qCurve;  //Total number of parameters
    
    fittedPara=gsl_matrix_alloc(cqsize, cnumOfPara);
    //To store the fitting result and error.
    fitErr=gsl_matrix_alloc(cqsize, cnumOfPara);
    status = new int[cqsize];		//Record the status of fitting.
    
    //Using Levenberg-Marquardt algorithm as implemented in the scaled lmder routine in minpack. Jacobian is given.
    const gsl_multifit_fdfsolver_type *solverType = gsl_multifit_fdfsolver_lmsder;
    
    int progress=0;		//Indicator of progress.
    
    //Objects to do numerical inverse Laplace transformation
#ifdef ISFRTD
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS), NILT3(OMP_NUM_THREADS), NILT4(OMP_NUM_THREADS);
#endif
    
#ifdef ISFRTDP
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS), NILT3(OMP_NUM_THREADS), NILT4(OMP_NUM_THREADS), NILT5(OMP_NUM_THREADS);
#endif
    
#ifdef ISFRTDPTT
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS), NILT3(OMP_NUM_THREADS), NILT4(OMP_NUM_THREADS), NILT5(OMP_NUM_THREADS), NILT6(OMP_NUM_THREADS);
#endif
    
#ifdef ISFRTDPfix
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS);
    
    const long double vbar=vbarGuess;
    const long double sigma=sigmaGuess;
    
    const long double vbsigma2=vbar/sigma/sigma;
    const long double vb2sigma2=vbsigma2*vbar;
    const long double logvbsigma2=log(vbsigma2);
    const long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    const long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    const long double vb2sigma3=vb2sigma2/sigma;
#endif
    
#ifdef ISFRTDPTTfix
    NILT NILT1(OMP_NUM_THREADS), NILT2(OMP_NUM_THREADS), NILT3(OMP_NUM_THREADS);
    
    const long double vbar=vbarGuess;
    const long double sigma=sigmaGuess;
    
    const long double vbsigma2=vbar/sigma/sigma;
    const long double vb2sigma2=vbsigma2*vbar;
    const long double logvbsigma2=log(vbsigma2);
    const long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    const long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    const long double vb2sigma3=vb2sigma2/sigma;
#endif
    
#pragma omp parallel for
    for (int iterq=0; iterq<cqsize; ++iterq)
    {
        //Data array which is going to present to the fitting algorithm
        double* datafit=new double[ctnum_fit];
        double* qList=new double[cnum_qCurve];
        double* time=new double[ctnum_fit];
        //Truncate the data, and put multiple curves into one array
        for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
        {
            for (int iterf = 0; iterf < cnum_fit; ++iterf)
            {
                datafit[iterf+iterqc*cnum_fit]=(gsl_matrix_get(datag, iterq+qIncreList[iterqc], iterf));		//Fitting in log scale.
                time[iterf+iterqc*cnum_fit]=tau[iterf];
            }
            qList[iterqc]=qabs[iterq+qIncreList[iterqc]];
        }
        
        gsl_multifit_function_fdf fitfun;		//Pointer of function to fit.
        dataStruct sdata;		//GSL data structure
        
        //Data is passed to ISFfun by sdata
        sdata.data=datafit;
        sdata.tau=time;
        sdata.q=qList;
        sdata.num_fit=cnum_fit;
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
        
#ifdef ISFRTDPTT
        sdata.ISFILT=&NILT1;
        sdata.dvbarISFILT=&NILT2;
        sdata.dsigmaISFILT=&NILT3;
        sdata.dDISFILT=&NILT4;
        sdata.dlambdaISFILT=&NILT5;
        sdata.dTTISFILT=&NILT6;
#endif
        
#ifdef ISFRTDPfix
        sdata.alpha=alphaGuess;
        sdata.D=DGuess;
        sdata.vbar=vbar;
        sdata.sigma=sigma;
        
        sdata.vbsigma2=vbsigma2;
        sdata.logfactor=logfactor;
        sdata.vb2sigma2=vb2sigma2;
        sdata.cpsiz1=cpsiz1;
        sdata.vb2sigma3=vb2sigma3;
        sdata.ISFILT=&NILT1;
        sdata.dlambdaISFILT=&NILT2;
#endif
        
#ifdef ISFRTDPTTfix
        sdata.alpha=alphaGuess;
        sdata.D=DGuess;
        sdata.vbar=vbar;
        sdata.sigma=sigma;
        
        sdata.vbsigma2=vbsigma2;
        sdata.logfactor=logfactor;
        sdata.vb2sigma2=vb2sigma2;
        sdata.cpsiz1=cpsiz1;
        sdata.vb2sigma3=vb2sigma3;
        sdata.ISFILT=&NILT1;
        sdata.dlambdaISFILT=&NILT2;
        sdata.dTTISFILT=&NILT3;
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
        fitfun.n=ctnum_fit;
        fitfun.p=cnumOfPara;
        fitfun.params=&sdata;
        
        //Initialization of the parameters
        double* localinipara=new double[cnumOfPara];
        for (int iterp=0; iterp<numOfPara; ++iterp)
        {
            localinipara[iterp]=inipara[iterp];
        }
        //Estimation of A(q) and B(q)
        for (int iterqc=0; iterqc<num_qCurve; ++iterqc)
        {
            localinipara[numOfPara+1+2*iterqc] = gsl_matrix_get(datag, iterq+qIncreList[iterqc], 0);
            localinipara[numOfPara+2*iterqc] = gsl_matrix_get(datag, iterq+qIncreList[iterqc], num_fit-1)-localinipara[numOfPara+1+2*iterqc];
        }
        //Initiallization of the solver
        gsl_vector_view para=gsl_vector_view_array(localinipara, cnumOfPara);
        gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(solverType, ctnum_fit, cnumOfPara);
        gsl_multifit_fdfsolver_set(solver, &fitfun, &para.vector);
        int iter=0;
        //gsl_vector* g=gsl_vector_alloc(numOfPara);
        
        //For debugging and monitering the iterations
//        cout << qList[0] << ' ' << qList[1] << '\n';
//        for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)
//        {
//            cout << gsl_vector_get(solver->x, iterpara) << '\n';
//        }
//        cout << '\n';
        
        int cstatus=GSL_CONTINUE;   //Current status
        do
        {
            gsl_multifit_fdfsolver_iterate(solver);		//Iterate one step.
            cstatus = norm0_rel_test(solver->dx, solver->x, 1e-7, 1e-7);  //Test the exiting criteria
            
            //For debugging and monitering the iterations
            //for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)
            //{
            //    cout << gsl_vector_get(solver->x, iterpara) << '\n';
            //}
            //cout << '\n';
            
            //If to use other exiting criteria
            //gsl_multifit_gradient(solver->J,solver->f, g);
            //status[iterq-1]=gsl_multifit_test_gradient(g, 1e-5);
            //			status[iterq - 1] = covar_rel_test(solver->J, solver->x, 1e-4);
            
            ++iter;
            //Number of iterations exceed certain limitation
            if (iter>maxIter)
            {
                cstatus=GSL_EMAXITER;
            }
        } while (cstatus == GSL_CONTINUE);
        status[iterq]=cstatus;
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
        
        //Output to standard I/O
        progress+=1;
        cout << "Fitted q=" << qabs[iterq] << " at iter=" << iter << ", " << 100.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << ", "<< gsl_strerror(status[iterq]) << "." << '\n';
        for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)
        {
            cout << gsl_matrix_get(fittedPara, iterq, iterpara) << '\n';
        }
        cout << '\n';
        delete [] datafit;
        delete [] qList;
        delete [] localinipara;
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