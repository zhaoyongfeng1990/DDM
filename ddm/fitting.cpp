//
//  fitting.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_multifit_nlin.h>
#include <omp.h>

int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J);

int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole);
//Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.

int covar_rel_test(const gsl_matrix* J, const gsl_vector* x, double tol);
//Test fitting result using error estimation from covariance matrix, not reliable. tol is the relative tolerance of error.

#ifndef MultiQFit
void ddm::fitting()
{
    fittedPara=gsl_matrix_alloc(qsize, numOfPara);
    //To store the fitting result and error.
    fitErr=gsl_matrix_alloc(qsize, numOfPara);
    status = new int[qsize];		//Record the status of fitting.
    datafit = gsl_matrix_alloc(qsize, num_fit);
    
    const gsl_multifit_fdfsolver_type *solverType;	//GSL solver
    solverType = gsl_multifit_fdfsolver_lmsder;
    //Using Levenberg-Marquardt algorithm as implemented in the scaled lmder routine in minpack. Jacobian is given.
    
    //Initial guess
    
    //Get the selected data for fitting
    for (int iterq = 0; iterq < qsize; ++iterq)
    {
        //int inif=floor(exp(-0.5)/dt);
        //int finalf=ceil(exp(0.5)/dt);
        for (int iterf = 0; iterf < num_fit; ++iterf)
        {
#ifdef NeedLaplaceTrans
            gsl_matrix_set(datafit, iterq, iterf, log(gsl_matrix_get(ldatag, iterq, iterf)));		//Fitting in log scale.
#else
            gsl_matrix_set(datafit, iterq, iterf, log(gsl_matrix_get(datag, iterq, iterf+iniTime)));		//Fitting in log scale.
#endif
        }
    }
    int progress=0;		//Indicator of progress.
    
#ifdef ISFRunAndTumbleAndDiffusionNoLT
    NILT NILT1, NILT2, NILT3, NILT4;
#endif
    
#ifdef ISFRTDPNoLT
    NILT NILT1, NILT2, NILT3, NILT4, NILT5;
#endif
    
#ifdef ISFRTDPNoLT_sigma
    NILT NILT1, NILT2, NILT3, NILT4, NILT5;
#endif
    
#pragma omp parallel for
    for (int iterq=0; iterq<qsize; ++iterq)
    {
        gsl_multifit_function_fdf fitfun;		//Function point.
        gsl_vector_view dataAry=gsl_matrix_row(datafit, iterq);
        dataStruct sdata;		//GSL data structure
        sdata.data=dataAry.vector.data;
#ifdef NeedLaplaceTrans
        sdata.tau=s;
#else
        sdata.tau=tau;
#endif
        sdata.q=qabs[iterq];
        sdata.num_fit=num_fit;
        
#ifdef ISFRunAndTumbleAndDiffusionNoLT
        sdata.ISFILT=&NILT1;
        sdata.dvISFILT=&NILT2;
        sdata.dDISFILT=&NILT3;
        sdata.dlambdaISFILT=&NILT4;
#endif
        
#ifdef ISFRTDPNoLT
        sdata.ISFILT=&NILT1;
        sdata.dvbarISFILT=&NILT2;
        sdata.dZISFILT=&NILT3;
        sdata.dDISFILT=&NILT4;
        sdata.dlambdaISFILT=&NILT5;
#endif
        
#ifdef ISFRTDPNoLT_sigma
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
        fitfun.n=num_fit;
        fitfun.p=numOfPara;
        fitfun.params=&sdata;
        
        double localinipara[numOfPara];
        for (int iterp=0; iterp<numOfPara-2; ++iterp)
        {
            localinipara[iterp]=inipara[iterp];
        }
        //Estimation of A(q) and B(q)
        localinipara[numOfPara-1] = gsl_matrix_get(datag, iterq, 0);
        localinipara[numOfPara-2] = gsl_matrix_get(datag, iterq, numOfDiff-1)-localinipara[numOfPara-1];
        
        //Initiallization of the solver
        gsl_vector_view para=gsl_vector_view_array(localinipara, numOfPara);
        gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(solverType, num_fit, numOfPara);
        gsl_multifit_fdfsolver_set(solver, &fitfun, &para.vector);
        int iter=0;
        //gsl_vector* g=gsl_vector_alloc(numOfPara);
        //        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        //        {
        //            cout << gsl_vector_get(solver->x, iterpara) << endl;
        //        }
        //        cout << endl;
        
        do
        {
            gsl_multifit_fdfsolver_iterate(solver);		//Iterate one step.
            status[iterq] = norm0_rel_test(solver->dx, solver->x, 1e-10, 1e-10);  //Test the exiting condition
            
            //            for (int iterpara=0; iterpara<numOfPara; ++iterpara)
            //            {
            //                cout << gsl_vector_get(solver->x, iterpara) << endl;
            //            }
            //            cout << endl;
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
        gsl_matrix* covar=gsl_matrix_alloc(numOfPara, numOfPara);
        gsl_multifit_covar(solver->J, 0.0, covar);
        for (int iterpara=0; iterpara<numOfPara; ++iterpara)	//Record result.
        {
            gsl_matrix_set(fittedPara, iterq, iterpara, gsl_vector_get(solver->x, iterpara) );
            gsl_matrix_set(fitErr, iterq, iterpara, sqrt(gsl_matrix_get(covar, iterpara, iterpara)) );    //Not presice in log scale
        }
        gsl_matrix_free(covar);
        gsl_multifit_fdfsolver_free(solver);
        
        progress+=1;
        cout << "Fitted q=" << qabs[iterq] << " at iter=" << iter << ", " << 100.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << ", "<< gsl_strerror(status[iterq]) << "." << endl;
        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        {
            cout << gsl_matrix_get(fittedPara, iterq, iterpara) << endl;
        }
    }
}

void ddm::fitting_estRange()
{
    fittedPara=gsl_matrix_alloc(qsize, numOfPara);
    //To store the fitting result and error.
    fitErr=gsl_matrix_alloc(qsize, numOfPara);
    status = new int[qsize];		//Record the status of fitting.
    //datafit = gsl_matrix_alloc(qsize, num_fit);
    
    const gsl_multifit_fdfsolver_type *solverType;	//GSL solver
    solverType = gsl_multifit_fdfsolver_lmsder;
    //Using Levenberg-Marquardt algorithm as implemented in the scaled lmder routine in minpack. Jacobian is given.
    
    int progress=0;		//Indicator of progress.
    
#ifdef ISFRunAndTumbleAndDiffusionNoLT
    NILT NILT1, NILT2, NILT3, NILT4;
#endif
    
#ifdef ISFRTDPNoLT
    NILT NILT1, NILT2, NILT3, NILT4, NILT5;
#endif
    
#ifdef ISFRTDPNoLT_sigma
    NILT NILT1, NILT2, NILT3, NILT4, NILT5;
#endif
    
#pragma omp parallel for
    for (int iterq=0; iterq<qsize; ++iterq)
    {
        double B = gsl_matrix_get(datag, iterq, 0);
        double A = gsl_matrix_get(datag, iterq, numOfDiff-1)-B;
        int ciniTime=-1;
        for (int itert=0; itert<numOfDiff; ++itert)
        {
            if (gsl_matrix_get(datag, iterq, itert)>B+0.5*A)
            {
                ciniTime=itert-1;
                break;
            }
        }
        if (ciniTime==-1)
        {
            progress+=1;
            cout << "Skipping q=" << qabs[iterq] << ", " << 100.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << "." << endl;
            continue;
        }
        int cfinalTime=10*ciniTime;
        cfinalTime=(cfinalTime>numOfDiff) ? numOfDiff : cfinalTime;
        ciniTime=0;
        int tempnum_fit=cfinalTime-ciniTime;
        double data[numOfDiff];
        double temptau[numOfDiff];
        for (int iterf = 0; iterf < tempnum_fit; ++iterf)
        {
            data[iterf]=log(gsl_matrix_get(datag, iterq, iterf+iniTime));		//Fitting in log scale.
            temptau[iterf]=(iterf+1+ciniTime)*dt;
        }
        
        gsl_multifit_function_fdf fitfun;		//Function point.
        dataStruct sdata;		//GSL data structure
        sdata.data=data;
        sdata.tau=temptau;
        sdata.q=qabs[iterq];
        sdata.num_fit=tempnum_fit;
        
#ifdef ISFRunAndTumbleAndDiffusionNoLT
        sdata.ISFILT=&NILT1;
        sdata.dvISFILT=&NILT2;
        sdata.dDISFILT=&NILT3;
        sdata.dlambdaISFILT=&NILT4;
#endif
        
#ifdef ISFRTDPNoLT
        sdata.ISFILT=&NILT1;
        sdata.dvbarISFILT=&NILT2;
        sdata.dZISFILT=&NILT3;
        sdata.dDISFILT=&NILT4;
        sdata.dlambdaISFILT=&NILT5;
#endif
        
#ifdef ISFRTDPNoLT_sigma
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
        fitfun.n=tempnum_fit;
        fitfun.p=numOfPara;
        fitfun.params=&sdata;
        
        //Estimation of A(q) and B(q)
        double localinipara[numOfPara];
        for (int iterp=0; iterp<numOfPara-2; ++iterp)
        {
            localinipara[iterp]=inipara[iterp];
        }
        //Estimation of A(q) and B(q)
        localinipara[numOfPara-1] = B;
        localinipara[numOfPara-2] = A;
        
        //Initiallization of the solver
        gsl_vector_view para=gsl_vector_view_array(localinipara, numOfPara);
        gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(solverType, tempnum_fit, numOfPara);
        gsl_multifit_fdfsolver_set(solver, &fitfun, &para.vector);
        int iter=0;
        //gsl_vector* g=gsl_vector_alloc(numOfPara);
        //        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        //        {
        //            cout << gsl_vector_get(solver->x, iterpara) << endl;
        //        }
        //        cout << endl;
        
        do
        {
            gsl_multifit_fdfsolver_iterate(solver);		//Iterate one step.
            status[iterq] = norm0_rel_test(solver->dx, solver->x, 1e-10, 1e-10);  //Test the exiting condition
            
            //            for (int iterpara=0; iterpara<numOfPara; ++iterpara)
            //            {
            //                cout << gsl_vector_get(solver->x, iterpara) << endl;
            //            }
            //            cout << endl;
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
        gsl_matrix* covar=gsl_matrix_alloc(numOfPara, numOfPara);
        gsl_multifit_covar(solver->J, 0.0, covar);
        for (int iterpara=0; iterpara<numOfPara; ++iterpara)	//Record result.
        {
            gsl_matrix_set(fittedPara, iterq, iterpara, gsl_vector_get(solver->x, iterpara) );
            gsl_matrix_set(fitErr, iterq, iterpara, sqrt(gsl_matrix_get(covar, iterpara, iterpara)) );    //Not presice in log scale
        }
        gsl_matrix_free(covar);
        gsl_multifit_fdfsolver_free(solver);
        
        progress+=1;
        cout << "Fitted q=" << qabs[iterq] << " at iter=" << iter << ", " << 100.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << ", "<< gsl_strerror(status[iterq]) << "." << endl;
        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        {
            cout << gsl_matrix_get(fittedPara, iterq, iterpara) << endl;
        }
    }
}
#endif

#ifdef MultiQFit
void ddm::fitting_DoubQ()
{
    int fqsize=floor(qsize/2);
    fittedPara=gsl_matrix_alloc(fqsize, numOfPara+2);
    //To store the fitting result and error.
    fitErr=gsl_matrix_alloc(fqsize, numOfPara+2);
    status = new int[fqsize];		//Record the status of fitting.
    //datafit = gsl_matrix_alloc(qsize, num_fit);
    double tempinipara[numOfPara+2];
    for (int iter=0; iter<numOfPara-2; ++iter)
    {
        tempinipara[iter]=inipara[iter];
    }
    
    const gsl_multifit_fdfsolver_type *solverType;	//GSL solver
    solverType = gsl_multifit_fdfsolver_lmsder;
    //Using Levenberg-Marquardt algorithm as implemented in the scaled lmder routine in minpack. Jacobian is given.
    
    int progress=0;		//Indicator of progress.
    
#ifdef ISFRunAndTumbleAndDiffusionNoLT
    NILT NILT1, NILT2, NILT3, NILT4;
#endif
    
#ifdef ISFRTDPNoLT
    NILT NILT1, NILT2, NILT3, NILT4, NILT5;
#endif
    
#ifdef ISFRTDPNoLT_sigma
    NILT NILT1, NILT2, NILT3, NILT4, NILT5;
#endif
    
#pragma omp parallel for
    for (int iterq=0; iterq<fqsize; ++iterq)
    {
        double B1 = gsl_matrix_get(datag, iterq, 0);
        double A1 = gsl_matrix_get(datag, iterq, numOfDiff-1)-B1;
        int iniTime1=-1;
        for (int itert=0; itert<numOfDiff; ++itert)
        {
            if (gsl_matrix_get(datag, iterq, itert)>B1+0.5*A1)
            {
                iniTime1=itert-1;
                break;
            }
        }
        if (iniTime1==-1)
        {
            progress+=1;
            cout << "Skipping q=" << qabs[iterq] << " and " << qabs[iterq*2] << ", " << 200.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << "." << endl;
            continue;
        }
        int finalTime1=10*iniTime1;
        finalTime1=(finalTime1>numOfDiff) ? numOfDiff : finalTime1;
        iniTime1=0;
        
        double B2 = gsl_matrix_get(datag, iterq*2, 0);
        double A2 = gsl_matrix_get(datag, iterq*2, numOfDiff-1)-B2;
        int iniTime2=-1;
        for (int itert=0; itert<numOfDiff; ++itert)
        {
            if (gsl_matrix_get(datag, iterq*2, itert)>B2+0.5*A2)
            {
                iniTime2=itert-1;
                break;
            }
        }
        if (iniTime2==-1)
        {
            progress+=1;
            cout << "Skipping q=" << qabs[iterq] << " and " << qabs[iterq*2] << ", " << 200.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << "." << endl;
            continue;
        }
        int finalTime2=10*iniTime2;
        finalTime2=(finalTime2>numOfDiff) ? numOfDiff : finalTime2;
        iniTime2=0;
        
        int num_fit1=finalTime1-iniTime1;
        int num_fit2=finalTime2-iniTime2;
        int tempnum_fit[2]={num_fit1, num_fit2};
        int num_fitt=num_fit1+num_fit2;
        
        double data[numOfDiff*2];
        double tempq[2]={qabs[iterq],qabs[iterq*2]};
        double tempt[numOfDiff*2];
        
        for (int iterf = 0; iterf < num_fit1; ++iterf)
        {
            data[iterf]=log(gsl_matrix_get(datag, iterq, iterf+iniTime1));		//Fitting in log scale.
            tempt[iterf]=(iterf+1+iniTime1)*dt;
        }
        for (int iterf = num_fit1; iterf < num_fitt; ++iterf)
        {
            data[iterf]=log(gsl_matrix_get(datag, iterq*2, iterf-num_fit1+iniTime2));		//Fitting in log scale.
            tempt[iterf]=(iterf-num_fit1+1+iniTime2)*dt;
        }
        
        gsl_multifit_function_fdf fitfun;		//Function point.
        dataStruct sdata;		//GSL data structure
        sdata.data=data;
        sdata.tau=tempt;
        sdata.q=tempq;
        sdata.num_fit=tempnum_fit;
        
#ifdef ISFRunAndTumbleAndDiffusionNoLT
        sdata.ISFILT=&NILT1;
        sdata.dvISFILT=&NILT2;
        sdata.dDISFILT=&NILT3;
        sdata.dlambdaISFILT=&NILT4;
#endif
        
#ifdef ISFRTDPNoLT
        sdata.ISFILT=&NILT1;
        sdata.dvbarISFILT=&NILT2;
        sdata.dZISFILT=&NILT3;
        sdata.dDISFILT=&NILT4;
        sdata.dlambdaISFILT=&NILT5;
#endif
        
#ifdef ISFRTDPNoLT_sigma
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
        fitfun.n=num_fitt;
        fitfun.p=numOfPara+2;
        fitfun.params=&sdata;
        
        //Estimation of A(q) and B(q)
        tempinipara[numOfPara-1] = B1;
        tempinipara[numOfPara-2] = A1;
        
        tempinipara[numOfPara+1] = B2;
        tempinipara[numOfPara] = A2;
        
        //Initiallization of the solver
        gsl_vector_view para=gsl_vector_view_array(tempinipara, numOfPara+2);
        gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(solverType, num_fitt, numOfPara+2);
        gsl_multifit_fdfsolver_set(solver, &fitfun, &para.vector);
        int iter=0;
        //gsl_vector* g=gsl_vector_alloc(numOfPara);
        //        for (int iterpara=0; iterpara<numOfPara; ++iterpara)
        //        {
        //            cout << gsl_vector_get(solver->x, iterpara) << endl;
        //        }
        //        cout << endl;
        
        do
        {
            gsl_multifit_fdfsolver_iterate(solver);		//Iterate one step.
            status[iterq] = norm0_rel_test(solver->dx, solver->x, 1e-10, 1e-10);  //Test the exiting condition
            
            //            for (int iterpara=0; iterpara<numOfPara; ++iterpara)
            //            {
            //                cout << gsl_vector_get(solver->x, iterpara) << endl;
            //            }
            //            cout << endl;
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
        gsl_matrix* covar=gsl_matrix_alloc(numOfPara+2, numOfPara+2);
        gsl_multifit_covar(solver->J, 0.0, covar);
        for (int iterpara=0; iterpara<numOfPara+2; ++iterpara)	//Record result.
        {
            gsl_matrix_set(fittedPara, iterq, iterpara, gsl_vector_get(solver->x, iterpara) );
            gsl_matrix_set(fitErr, iterq, iterpara, sqrt(gsl_matrix_get(covar, iterpara, iterpara)) );    //Not presice in log scale
        }
        gsl_matrix_free(covar);
        gsl_multifit_fdfsolver_free(solver);
        
        progress+=1;
        cout << "Fitted q=" << qabs[iterq] << " at iter=" << iter << ", " << 100.0*progress / qsize << "% completed from thread No." << omp_get_thread_num() << ", "<< gsl_strerror(status[iterq]) << "." << endl;
        for (int iterpara=0; iterpara<numOfPara+2; ++iterpara)
        {
            cout << gsl_matrix_get(fittedPara, iterq, iterpara) << endl;
        }
    }
}
#endif

#ifndef NoJacobian
int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J)
{
    ISFfun(para, sdata, y);
    dISFfun(para, sdata, J);
    
    return  GSL_SUCCESS;
}
#endif

//Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.
#ifndef MultiQFit
int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole)
{
    //Check if some parameters become NAN or INF.
    for (int iter = 0; iter < numOfPara; ++iter)
    {
        if (!isfinite(x->data[iter]))
        {
            return GSL_EOVRFLW;
        }
    }
    for (int iter = 0; iter < numOfPara; ++iter)
    {
        //double relerr = abs(dx->data[iter]/* / x->data[iter]*/);
        if (abs(dx->data[iter]) > tol*abs(x->data[iter])+tole)
        {
            return GSL_CONTINUE;
        }
    }
    return GSL_SUCCESS;
}
#else
int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole)
{
    //Check if some parameters become NAN or INF.
    for (int iter = 0; iter < numOfPara+2; ++iter)
    {
        if (!isfinite(x->data[iter]))
        {
            return GSL_EOVRFLW;
        }
    }
    for (int iter = 0; iter < numOfPara+2; ++iter)
    {
        //double relerr = abs(dx->data[iter]/* / x->data[iter]*/);
        if (abs(dx->data[iter]) > tol*abs(x->data[iter])+tole)
        {
            return GSL_CONTINUE;
        }
    }
    return GSL_SUCCESS;
}
#endif

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