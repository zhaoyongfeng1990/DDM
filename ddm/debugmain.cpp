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
    long double q=3.00l;
    long double v0=20.0l;//gsl_vector_get(para, 1);
    long double lambda=100.0l;
    long double D=0.2l;//gsl_vector_get(para, 3);\

    long double kv0=q*v0;
    long double Dq2=D*q*q;
    long double Dq2lambda=Dq2+lambda;
    long double qv2=kv0*kv0;
    long double paraISF[6]={kv0, lambda, Dq2, q, Dq2lambda, qv2};
    
    int tid=omp_get_thread_num();
    
    long double qvlambda=kv0/lambda;
    long double incre=0.01l;
    
    if (qvlambda>(pi/2))
    {
        ILT.weideman(-Dq2lambda, kv0, incre);
        dvILT.weideman(-Dq2lambda, kv0, incre);
        dlambdaILT.weideman(-Dq2lambda, kv0, incre);
        dDILT.weideman(-Dq2lambda, kv0, incre);
    }
    else
    {
        ILT.weideman(-Dq2lambda, kv0, -Dq2lambda+kv0/tan(qvlambda),0, incre);
        dvILT.weideman(-Dq2lambda, kv0, -Dq2lambda+kv0/tan(qvlambda),0, incre);
        dlambdaILT.weideman(-Dq2lambda, -Dq2lambda+kv0/tan(qvlambda),0, kv0, incre);
        dDILT.weideman(-Dq2lambda, kv0, -Dq2lambda+kv0/tan(qvlambda),0, incre);
    }
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