//
//  debugmain.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/5.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include "NILT.h"
#include "ddm.h"
using namespace std;
//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_sf_psi.h>
//
//gsl_integration_cquad_workspace* workspace;
//
cpx ISFs(cpx s, long double* para, long double v);
cpx dvbarISFs(cpx s, long double* para, long double v);
cpx dsigmaISFs(cpx s, long double* para, long double v);
cpx dDISFs(cpx s, long double* para, long double v);
cpx dlambdaISFs(cpx s, long double* para, long double v);

int main()
{
    NILT ILT(1);
    NILT dvbarILT(1);
    NILT dsigmaILT(1);
    NILT dlambdaILT(1);
    NILT dDILT(1);
    long double q=0.48143l;
    long double vbar=12.569l;//gsl_vector_get(para, 1);
    long double sigma=4.2482;
    long double lambda=3.3242e-8l;
    long double D=0.55492l;//gsl_vector_get(para, 3);\
   
    const long double Dq2=D*q*q;
    const long double Dq2lambda=Dq2+lambda;
    const long double vbsigma2=vbar/sigma/sigma;
    const long double vb2sigma2=vbsigma2*vbar;
    const long double logvbsigma2=log(vbsigma2);
    const long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    const long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    const long double vb2sigma3=vb2sigma2/sigma;
    
    long double paraISF[10]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2, vbar, cpsiz1, vb2sigma3, sigma};
    
    int tid=0;
    
    ILT.cfun[tid].fun=ISFs;
    dvbarILT.cfun[tid].fun=dvbarISFs;
    dsigmaILT.cfun[tid].fun=dsigmaISFs;
    dlambdaILT.cfun[tid].fun=dlambdaISFs;
    dDILT.cfun[tid].fun=dDISFs;
    
    long double& csigma=ILT.sigma[tid];
    long double& cb=ILT.b[tid];
    long double& cb2=ILT.b2[tid];
    long double& csigmab=ILT.sigmab[tid];
    
    csigma=-Dq2+0.5l;
    cb=sqrt(csigma*csigma-Dq2*Dq2+(Dq2+Dq2lambda)*(csigma+Dq2));
    cb2=cb*2;
    csigmab=csigma-cb;
    
    dvbarILT.sigma[tid]=csigma;
    dsigmaILT.sigma[tid]=csigma;
    dlambdaILT.sigma[tid]=csigma;
    dDILT.sigma[tid]=csigma;
    
    dvbarILT.b[tid]=cb;
    dsigmaILT.b[tid]=cb;
    dlambdaILT.b[tid]=cb;
    dDILT.b[tid]=cb;
    
    dvbarILT.b2[tid]=cb2;
    dsigmaILT.b2[tid]=cb2;
    dlambdaILT.b2[tid]=cb2;
    dDILT.b2[tid]=cb2;
    
    dvbarILT.sigmab[tid]=csigmab;
    dsigmaILT.sigmab[tid]=csigmab;
    dlambdaILT.sigmab[tid]=csigmab;
    dDILT.sigmab[tid]=csigmab;
    
    ILT.NiLT_weeks(paraISF);
    dvbarILT.NiLT_weeks(paraISF);
    dsigmaILT.NiLT_weeks(paraISF);
    dDILT.NiLT_weeks(paraISF);
    dlambdaILT.NiLT_weeks(paraISF);
    
    //    cb2=cb*2;
    //    csigmab=csigma-cb;
    
    ofstream debug("debug.txt");
    ofstream dvdebug("dvdebug.txt");
    ofstream dZdebug("dZdebug.txt");
    ofstream dlambdadebug("dlambdadebug.txt");
    ofstream dDdebug("dDdebug.txt");
    //ofstream a("a.txt");
    //cout << M << '\n';
    //for (int iter = 0; iter<M; ++iter)
    //{
    //a << setprecision(30) << ILT.CoeA[tid][iter] << '\n';
    //Temperary variables used for acceleration.
    //}
    for (int iter = 0; iter<9000; ++iter)
    {
        long double t=(iter+1)*0.01l;
        double rtd=ILT.clenshaw(t);
        double dvrtd=dvbarILT.clenshaw(t);
        double dsigmartd=dsigmaILT.clenshaw(t);
        double dlambdartd=dlambdaILT.clenshaw(t);
        double dDrtd=dDILT.clenshaw(t);
        debug << setprecision(30) << rtd << '\n';
        dvdebug << setprecision(30) << dvrtd << '\n';
        dlambdadebug << setprecision(30) << dlambdartd << '\n';
        dDdebug << setprecision(30) << dDrtd << '\n';
        dZdebug << setprecision(30) << dsigmartd << '\n';
        //Temperary variables used for acceleration.
    }
    
    return 0;
}