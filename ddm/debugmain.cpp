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
#include "NILT.h"
#include "ddm.h"
using namespace std;
//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_sf_psi.h>
//
//gsl_integration_cquad_workspace* workspace;
//
cpx ISFs(cpx s, long double* para);
cpx dvISFs(cpx s, long double* para);
cpx dDISFs(cpx s, long double* para);
cpx dlambdaISFs(cpx s, long double* para);

int main()
{
    NILT ILT;
    NILT dvILT;
    NILT dlambdaILT;
    NILT dDILT;
    long double q=0.02l;
    long double v0=11.2917l;//gsl_vector_get(para, 1);
    long double lambda=0.829098l;
    long double D=9.0517l;//gsl_vector_get(para, 3);\

    long double kv0=q*v0;
    long double Dq2=D*q*q;
    long double Dq2lambda=Dq2+lambda;
    long double qv2=kv0*kv0;
    long double paraISF[6]={kv0, lambda, Dq2, q, Dq2lambda, qv2};
    
    long double qvlambda=kv0/lambda;
    long double incre=0.01l;
    
    int tid=omp_get_thread_num();
    long double& csigma=ILT.sigma[tid];
    long double& cb=ILT.b[tid];
    long double& cb2=ILT.b2[tid];
    long double& csigmab=ILT.sigmab[tid];
    
    if (qvlambda>(pi/2))
    {
        csigma=-Dq2lambda+1.0l;
        cb=sqrt(kv0*kv0+1.0l*1.0l);
    }
    else
    {
        long double alpha21=kv0/tan(qvlambda);
        long double alpha1=-Dq2lambda;
        long double alpha2=alpha21+alpha1;
        csigma=alpha2+0.001l;
        cb=sqrt(csigma*csigma-((csigma-alpha1)*alpha2*alpha2-0.001l*(alpha1*alpha1+kv0*kv0))/alpha21);
    }
    cb2=cb*2;
    csigmab=csigma-cb;
    
    dvILT.sigma[tid]=csigma;
    dlambdaILT.sigma[tid]=csigma;
    dDILT.sigma[tid]=csigma;
    
    dvILT.b[tid]=cb;
    dlambdaILT.b[tid]=cb;
    dDILT.b[tid]=cb;
    
    dvILT.b2[tid]=cb2;
    dlambdaILT.b2[tid]=cb2;
    dDILT.b2[tid]=cb2;
    
    dvILT.sigmab[tid]=csigmab;
    dlambdaILT.sigmab[tid]=csigmab;
    dDILT.sigmab[tid]=csigmab;

//    cb2=cb*2;
//    csigmab=csigma-cb;
    
    ILT.NiLT_weeks(ISFs, paraISF);
    dvILT.NiLT_weeks(dvISFs, paraISF);
    dlambdaILT.NiLT_weeks(dlambdaISFs, paraISF);
    dDILT.NiLT_weeks(dDISFs, paraISF);
    ofstream debug("debug.txt");
    ofstream dvdebug("dvdebug.txt");
    ofstream dlambdadebug("dlambdadebug.txt");
    ofstream dDdebug("dDdebug.txt");
    ofstream a("a.txt");
    cout << M << endl;
    for (int iter = 0; iter<M; ++iter)
    {
      a << setprecision(30) << ILT.CoeA[tid][iter] << endl;
      //Temperary variables used for acceleration.
    }
    for (int iter = 0; iter<num_fit; ++iter)
    {
        long double t=(iter+1)*0.01l;
        double rtd=ILT.clenshaw(t);
        double dvrtd=dvILT.clenshaw(t);
        double dlambdartd=dlambdaILT.clenshaw(t);
        double dDrtd=dDILT.clenshaw(t);
        debug << setprecision(30) << rtd << endl;
        dvdebug << setprecision(30) << dvrtd << endl;
        dlambdadebug << setprecision(30) << dlambdartd << endl;
        dDdebug << setprecision(30) << dDrtd << endl;
        //Temperary variables used for acceleration.
    }

    return 0;
}