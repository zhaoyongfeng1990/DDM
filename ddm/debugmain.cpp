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
cpx dZISFs(cpx s, long double* para, long double v);
cpx dDISFs(cpx s, long double* para, long double v);
cpx dlambdaISFs(cpx s, long double* para, long double v);

int main()
{
    NILT ILT;
    NILT dvbarILT;
    NILT dZILT;
    NILT dlambdaILT;
    NILT dDILT;
    long double q=0.52l;
    long double vbar=13.0541l;//gsl_vector_get(para, 1);
    long double Z=127;
    long double lambda=1.07861l;
    long double D=0.399742l;//gsl_vector_get(para, 3);\
    
    long double Dq2=D*q*q;
    long double Dq2lambda=Dq2+lambda;
    long double Z1=Z+1.0l;
    long double Z1vbar=Z1/vbar;
    long double logfactor=Z1*log(Z1vbar)-gsl_sf_lngamma(Z1);
    long double psiz1=gsl_sf_psi(Z1);
    
    long double paraISF[8]={lambda, q, Dq2lambda, Z, logfactor, Z1vbar, vbar, psiz1};
    
    
    int tid=omp_get_thread_num();
    
    ILT.cfun[tid].fun=ISFs;
    dvbarILT.cfun[tid].fun=dvbarISFs;
    dZILT.cfun[tid].fun=dZISFs;
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
    dZILT.sigma[tid]=csigma;
    dlambdaILT.sigma[tid]=csigma;
    dDILT.sigma[tid]=csigma;
    
    dvbarILT.b[tid]=cb;
    dZILT.b[tid]=cb;
    dlambdaILT.b[tid]=cb;
    dDILT.b[tid]=cb;
    
    dvbarILT.b2[tid]=cb2;
    dZILT.b2[tid]=cb2;
    dlambdaILT.b2[tid]=cb2;
    dDILT.b2[tid]=cb2;
    
    dvbarILT.sigmab[tid]=csigmab;
    dZILT.sigmab[tid]=csigmab;
    dlambdaILT.sigmab[tid]=csigmab;
    dDILT.sigmab[tid]=csigmab;
    
    //ILT.NiLT_weeks(paraISF);
    //dvbarILT.NiLT_weeks(paraISF);
    //dZILT.NiLT_weeks(paraISF);
    //dDILT.NiLT_weeks(paraISF);
    dlambdaILT.NiLT_weeks(paraISF);

//    cb2=cb*2;
//    csigmab=csigma-cb;
    
    ofstream debug("debug.txt");
    ofstream dvdebug("dvdebug.txt");
    ofstream dZdebug("dZdebug.txt");
    ofstream dlambdadebug("dlambdadebug.txt");
    ofstream dDdebug("dDdebug.txt");
    ofstream a("a.txt");
    cout << M << '\n';
    for (int iter = 0; iter<M; ++iter)
    {
      a << setprecision(30) << ILT.CoeA[tid][iter] << '\n';
      //Temperary variables used for acceleration.
    }
    for (int iter = 0; iter<num_fit; ++iter)
    {
        long double t=(iter+1)*0.01l;
        double rtd=ILT.clenshaw(t);
        double dvrtd=dvbarILT.clenshaw(t);
        double dZrtd=dZILT.clenshaw(t);
        double dlambdartd=dlambdaILT.clenshaw(t);
        double dDrtd=dDILT.clenshaw(t);
        debug << setprecision(30) << rtd << '\n';
        dvdebug << setprecision(30) << dvrtd << '\n';
        dlambdadebug << setprecision(30) << dlambdartd << '\n';
        dDdebug << setprecision(30) << dDrtd << '\n';
        dZdebug << setprecision(30) << dZrtd << '\n';
        //Temperary variables used for acceleration.
    }

    return 0;
}