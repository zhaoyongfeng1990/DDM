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
#include <gsl/gsl_sf_lambert.h>

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
cpx dTTISFs(cpx s, long double* para, long double v);

int main()
{
    NILT ILT(1);
    NILT dvbarILT(1);
    NILT dsigmaILT(1);
    NILT dlambdaILT(1);
    NILT dDILT(1);
    NILT dTTILT(1);
    long double q=0.2l;
    long double vbar=30.0l;//gsl_vector_get(para, 1);
    long double sigma=6.1l;
    long double lambda=0.001l;
    long double D=0.1l;//gsl_vector_get(para, 3);
    long double TT=1e-3l;
    
    const long double Dq2=D*q*q;
    const long double Dq2lambda=Dq2+lambda;
    const long double vbsigma2=vbar/sigma/sigma;
    const long double vb2sigma2=vbsigma2*vbar;
    const long double logvbsigma2=log(vbsigma2);
    const long double logfactor=vb2sigma2*logvbsigma2-gsl_sf_lngamma(vb2sigma2);
    const long double cpsiz1=logvbsigma2-gsl_sf_psi(vb2sigma2);
    const long double vb2sigma3=vb2sigma2/sigma;
    
    const long double tumbFrac=1.0l+lambda*TT;
    const long double lambdaTT=lambda*TT;
    const long double lambdaTT2=lambdaTT*TT;
    const long double TT2=TT*TT;
    const long double lambda2=lambda*lambda;
    const long double lambda3=lambda2*lambda;
    
    long double paraISF[18]={lambda, q, Dq2lambda, vbsigma2, logfactor, vb2sigma2, vbar, cpsiz1, vb2sigma3, sigma, TT, Dq2, tumbFrac, lambdaTT2, lambdaTT, TT2, lambda2, lambda3};
    
    int tid=0;
    
    ILT.cfun[tid].fun=ISFs;
    dvbarILT.cfun[tid].fun=dvbarISFs;
    dsigmaILT.cfun[tid].fun=dsigmaISFs;
    dlambdaILT.cfun[tid].fun=dlambdaISFs;
    dDILT.cfun[tid].fun=dDISFs;
    dTTILT.cfun[tid].fun=dTTISFs;
    
    long double& csigma=ILT.sigma[tid];
    long double& cb=ILT.b[tid];
    long double& cb2=ILT.b2[tid];
    long double& csigmab=ILT.sigmab[tid];
    
    long double a1=-Dq2-lambda-1/TT;
    long double a21=lambda+1/TT;
    
    long double a2=-Dq2;
    
    csigma=a2+0.5l;
    
    cb=sqrt(csigma*csigma-(a2*a2*(csigma-a1)-a1*a1*(csigma-a2))/a21);
    
    //csigma=-Dq2+0.5l;
    //cb=sqrt(csigma*csigma-(csigma*(1.0l+lambda*TT)+(Dq2*TT/(1.0l+lambda*TT)+1.0l)*Dq2*(1.0l+lambda*TT+2.0l*sigma*TT))/TT);
    cb2=cb*2;
    csigmab=csigma-cb;
    
    dvbarILT.sigma[tid]=csigma;
    dsigmaILT.sigma[tid]=csigma;
    dlambdaILT.sigma[tid]=csigma;
    dDILT.sigma[tid]=csigma;
    dTTILT.sigma[tid]=csigma;
    
    dvbarILT.b[tid]=cb;
    dsigmaILT.b[tid]=cb;
    dlambdaILT.b[tid]=cb;
    dDILT.b[tid]=cb;
    dTTILT.b[tid]=cb;
    
    dvbarILT.b2[tid]=cb2;
    dsigmaILT.b2[tid]=cb2;
    dlambdaILT.b2[tid]=cb2;
    dDILT.b2[tid]=cb2;
    dTTILT.b2[tid]=cb2;
    
    dvbarILT.sigmab[tid]=csigmab;
    dsigmaILT.sigmab[tid]=csigmab;
    dlambdaILT.sigmab[tid]=csigmab;
    dDILT.sigmab[tid]=csigmab;
    dTTILT.sigmab[tid]=csigmab;
    
    cout << "begin" << '\n';
    ILT.NiLT_weeks(paraISF);
    cout << "ISF ok" << '\n';
    dvbarILT.NiLT_weeks(paraISF);
    cout << "dvISF ok" << '\n';
    dsigmaILT.NiLT_weeks(paraISF);
    cout << "dsISF ok" << '\n';
    dDILT.NiLT_weeks(paraISF);
    cout << "ddISF ok" << '\n';
    dlambdaILT.NiLT_weeks(paraISF);
    cout << "dlISF ok" << '\n';
    dTTILT.NiLT_weeks(paraISF);
    cout << "dtISF ok" << '\n';
    
    //    cb2=cb*2;
    //    csigmab=csigma-cb;
    
    ofstream debug("debug.txt");
    ofstream dvdebug("dvdebug.txt");
    ofstream dZdebug("dZdebug.txt");
    ofstream dlambdadebug("dlambdadebug.txt");
    ofstream dDdebug("dDdebug.txt");
    ofstream dTTdebug("dTTdebug.txt");
    ofstream a("a.txt");
    cout << M << '\n';
    for (int iter = 0; iter<M; ++iter)
    {
        a << setprecision(30) << ILT.CoeA[tid][iter] << '\n';
    //Temperary variables used for acceleration.
    }
    for (int iter = 0; iter<4000; ++iter)
    {
        long double t=(iter+1)*0.01l;
        double rtd=ILT.clenshaw(t);
        double dvrtd=dvbarILT.clenshaw(t);
        double dsigmartd=dsigmaILT.clenshaw(t);
        double dlambdartd=dlambdaILT.clenshaw(t);
        double dDrtd=dDILT.clenshaw(t);
        double dTTrtd=dTTILT.clenshaw(t);
        debug << setprecision(30) << rtd << '\n';
        dvdebug << setprecision(30) << dvrtd << '\n';
        dlambdadebug << setprecision(30) << dlambdartd << '\n';
        dDdebug << setprecision(30) << dDrtd << '\n';
        dZdebug << setprecision(30) << dsigmartd << '\n';
        dTTdebug << setprecision(30) << dTTrtd << '\n';
        //Temperary variables used for acceleration.
    }
    
    return 0;
}