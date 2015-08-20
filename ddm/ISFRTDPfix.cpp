//
//  ISFRTDPfix.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/7/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
//#include <iostream>

#ifdef ISFRTDPfix

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <vector>
#include <omp.h>

//cpx is the complex long double datatype

//The Laplace transformed ISFs and its derivatives.
cpx ISFs(cpx s, long double* para, long double v);
cpx dlambdaISFs(cpx s, long double* para, long double v);

cpx ISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    const long double lambda=para[0];
    const long double q=para[1];
    const long double Dq2lambda=para[2];
    const long double vbsigma2=para[3];
    const long double logfactor=para[4];
    const long double vb2sigma2=para[5];
    
    const long double qv=q*v;
    const long double pv=exp((vb2sigma2-1)*log(v)+logfactor-v*vbsigma2);
    
    const cpx atanterm=atan(qv/(Dq2lambda+s));
    
    return atanterm/(qv-lambda*atanterm)*pv;
}

cpx dlambdaISFs(cpx s, long double* para, long double v)
{
    //Temperary variables used for acceleration.
    const long double lambda=para[0];
    const long double q=para[1];
    const long double Dq2lambda=para[2];
    const long double vbsigma2=para[3];
    const long double logfactor=para[4];
    const long double vb2sigma2=para[5];
    
    const long double qv=q*v;
    const long double pv=exp((vb2sigma2-1)*log(v)+logfactor-v*vbsigma2);
    
    const long double qv2=qv*qv;
    
    const cpx Dq2lambdas=Dq2lambda+s;
    const cpx atanterm=atan(qv/Dq2lambdas);
    return ((qv2+Dq2lambdas*Dq2lambdas)*atanterm*atanterm-qv2)/(qv2+Dq2lambdas*Dq2lambdas)/(qv-lambda*atanterm)/(qv-lambda*atanterm)*pv;
}

//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    //data list
    const double* dataAry=((dataStruct *)sdata)->data;
    //Time points list
    const double* tau=((dataStruct *)sdata)->tau;
    //q list
    const double* qArray=((dataStruct *)sdata)->q;
    //Get the class for numerical inverse Laplace Transformation
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    
    //Get the number of data points of each curve
    const int num_fit=((dataStruct *)sdata)->num_fit;
    //Get the number of curves
    const int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    //Get the total number of data points
    const int tnum_fit=num_fit*cnum_qCurve;
    //Get the parameters.
    const long double lambda=gsl_vector_get(para, 0);
    
    const long double alpha=((dataStruct *)sdata)->alpha;
    const long double D=((dataStruct *)sdata)->D;
    const long double vbar=((dataStruct *)sdata)->vbar;
    const long double sigma=((dataStruct *)sdata)->sigma;
    
    const long double vbsigma2=((dataStruct *)sdata)->vbsigma2;
    const long double logfactor=((dataStruct *)sdata)->logfactor;
    const long double vb2sigma2=((dataStruct *)sdata)->vb2sigma2;
    
    //    cout << qArray[0] << ' ' << qArray[1] << '\n';
    //    for (int i=0; i<9; ++i)
    //    {
    //        cout << gsl_vector_get(para, i) << '\n';
    //    }
    //    cout << '\n';
    
    //If the parameters are out of range, this flag will be turned on and a punishment term will be given as return.
    bool breakFlag=false;
    
    //Punishment terms, to make constrains in parameter space.
    double punishment=0;
    if (lambda<0)
    {
        punishment+=1e5*lambda*lambda;
        breakFlag=true;
    }
    if (lambda>1e6)
    {
        punishment+=1e5*(lambda-1e6)*(lambda-1e6);
        breakFlag=true;
    }
    
    if (breakFlag)
    {
        for (int iter = 0; iter<tnum_fit; ++iter)
        {
            gsl_vector_set(y, iter, punishment - dataAry[iter]);
        }
        return GSL_SUCCESS;
    }
    
    //Temperary variables used for acceleration.
    
    long double paraISF[10]={lambda, 0, 0, vbsigma2, logfactor, vb2sigma2, vbar, 0, 0, sigma};
    
    //The reference for setting the parameters of iLT solver
    int tid=omp_get_thread_num();
    
    ILT->cfun[tid].fun=ISFs;
    
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    //Loop over the curves
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        //Each curve have different value of q, but A(q) and B(q) depend on q.
        const long double A=gsl_vector_get(para, 1+2*iterqc);
        const long double B=gsl_vector_get(para, 2+2*iterqc);
        
        //Check if A is out of range
        if (A<0)
        {
            punishment+=1e5*A*A;
            breakFlag=true;
        }
        if (breakFlag)
        {
            for (int iter = 0; iter<tnum_fit; ++iter)
            {
                gsl_vector_set(y, iter, punishment - dataAry[iter]);
            }
            return GSL_SUCCESS;
        }
        
        //Temperary variables used for acceleration. I don't pass the parameters directly to the solver, because some terms need only be calculated only once.
        const long double q=qArray[iterqc];
        const long double Dq2=D*q*q;
        const double Dq2lambda=Dq2+lambda;
        
        paraISF[1]=q;
        paraISF[2]=Dq2lambda;
        
        //Initialization of numerical inverse Laplace transformation solver
        //Set sigma and b in iLT solver, using weideman's method. (To speed up, I didn't call function weidman(). )
        csigma=-Dq2+0.5l;
        cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
        
        //Temperary variables used for acceleration.
        cb2=cb*2;
        csigmab=csigma-cb;
        
        //Calculate the coefficients of Laguerre polynomial series expansion.
        ILT->NiLT_weeks(paraISF);
        
        //Loop over each data point of the curve
        for (int iter = 0; iter<num_fit; ++iter)
        {
            //The real index of iter-th data in each curve
            const int cidx=iter+iterqc*num_fit;
            const long double t=tau[cidx];
            //Evaluate ISF at time t, the coefficients has been calculated.
            const double rtd=ILT->clenshaw(t);
            const double yi=(A*(1.0-(1.0-alpha)*exp(-Dq2*t)-alpha*rtd)+B);
            
            //Actually, sqrt(weight)
            const double weight=1.0/sqrt(dataAry[cidx]);
            double result = (yi - dataAry[cidx])*weight;
            gsl_vector_set(y, cidx, result);
        }
    }
    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    //Get the number of data points of each curve
    const int num_fit=((dataStruct *)sdata)->num_fit;
    //Get the number of curves
    const int cnum_qCurve=((dataStruct *)sdata)->num_qCurve;
    //Get the total number of data points
    const int tnum_fit=num_fit*cnum_qCurve;
    
    //Get the parameters.
    const long double lambda=gsl_vector_get(para, 0);
    
    const long double alpha=((dataStruct *)sdata)->alpha;
    const long double D=((dataStruct *)sdata)->D;
    const long double vbar=((dataStruct *)sdata)->vbar;
    const long double sigma=((dataStruct *)sdata)->sigma;
    
    const long double vbsigma2=((dataStruct *)sdata)->vbsigma2;
    const long double logfactor=((dataStruct *)sdata)->logfactor;
    const long double vb2sigma2=((dataStruct *)sdata)->vb2sigma2;
    const long double cpsiz1=((dataStruct *)sdata)->cpsiz1;
    const long double vb2sigma3=((dataStruct *)sdata)->vb2sigma3;
    
    //Cleaning
    gsl_matrix_set_zero(J);
    //If the parameters are out of range, this flag will be turned on and a punishment term will be given as return.
    bool breakFlag=false;
    //Punishment terms, to make constrains in parameter space.
    if (lambda<0)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 2e5*lambda);
        }
        breakFlag=true;
    }
    if (lambda>1e6)
    {
        for (int iter=0; iter<tnum_fit; ++iter)
        {
            gsl_matrix_set(J, iter, 0, 2e5*(lambda-1e6));
        }
        breakFlag=true;
    }
    if (breakFlag)
    {
        return GSL_SUCCESS;
    }
    
    //Time points list
    const double* tau=((dataStruct *)sdata)->tau;
    //q list
    const double* qArray=((dataStruct *)sdata)->q;
    //data list
    const double* dataAry=((dataStruct *)sdata)->data;
    
    //Temperary variables used for acceleration.
    
    long double paraISF[10]={lambda, 0, 0, vbsigma2, logfactor, vb2sigma2, vbar, cpsiz1, vb2sigma3, sigma};
    
    //Get the class for numerical inverse Laplace Transformation, one for each function
    NILT* ILT=((dataStruct *)sdata)->ISFILT;
    NILT* dlambdaILT=((dataStruct *)sdata)->dlambdaISFILT;
    
    const int tid=omp_get_thread_num();
    
    ILT->cfun[tid].fun=ISFs;
    dlambdaILT->cfun[tid].fun=dlambdaISFs;
    
    //The reference for setting the parameters of iLT solver
    long double& csigma=ILT->sigma[tid];
    long double& cb=ILT->b[tid];
    long double& cb2=ILT->b2[tid];
    long double& csigmab=ILT->sigmab[tid];
    
    //Loop over the curves
    for (int iterqc=0; iterqc<cnum_qCurve; ++iterqc)
    {
        //Each curve have different value of q, but A(q) and B(q) depend on q.
        const long double A=gsl_vector_get(para, 1+2*iterqc);
        const long double B=gsl_vector_get(para, 2+2*iterqc);
        
        //Check if A is out of range
        if (A<0)
        {
            for (int iter=0; iter<tnum_fit; ++iter)
            {
                gsl_matrix_set(J, iter, 1+2*iterqc, 2e5*A);
            }
            breakFlag=true;
        }
        
        if (breakFlag)
        {
            return GSL_SUCCESS;
        }
        
        //Temperary variables used for acceleration. I don't pass the parameters directly to the solver, because some terms need only be calculated only once.
        const long double q=qArray[iterqc];
        const long double Dq2=D*q*q;
        const long double Dq2lambda=Dq2+lambda;
        
        paraISF[1]=q;
        paraISF[2]=Dq2lambda;
        
        //Initialization of numerical inverse Laplace transformation solver
        //Set sigma and b in iLT solver, using weideman's method. (To speed up, I didn't call function weidman(). )
        csigma=-Dq2+0.5l;
        cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
        
        cb2=cb*2;
        csigmab=csigma-cb;
        
        //All the functions have the same singularities, so all the solvers have the same parameters
        dlambdaILT->sigma[tid]=csigma;
        dlambdaILT->b[tid]=cb;
        dlambdaILT->b2[tid]=cb2;
        dlambdaILT->sigmab[tid]=csigmab;
        
        //Calculate the coefficients of Laguerre polynomial series expansion.
        ILT->NiLT_weeks(paraISF);
        dlambdaILT->NiLT_weeks(paraISF);
        
        //Loop over each data point of the curve
        for (int iter=0; iter<num_fit; ++iter)
        {
            //The real index of iter-th data in each curve
            const int cidx=iter+iterqc*num_fit;
            //Temperary variables used for acceleration.
            const long double t=tau[cidx];
            //Evaluate ISF at time t, the coefficients has been calculated.
            const double rtd=ILT->clenshaw(t);
            const double dlambdartd=dlambdaILT->clenshaw(t);
            
            const double expterm=exp(-Dq2*t);
            const double dA=(1.0-(1.0-alpha)*expterm-alpha*rtd);
            //Actually, sqrt(weight)
            const double weight=1.0/sqrt(dataAry[cidx]);
            
            gsl_matrix_set(J, cidx, 0, -A*alpha*dlambdartd*weight );
            gsl_matrix_set(J, cidx, 1+2*iterqc, dA*weight );
            gsl_matrix_set(J, cidx, 2+2*iterqc, 1.0*weight );
        }
    }
    return GSL_SUCCESS;
}

#endif
