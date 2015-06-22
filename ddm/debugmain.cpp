////
////  debugmain.cpp
////  ddm
////
////  Created by Zhao Yongfeng on 15/6/5.
////  Copyright (c) 2015年 ZYF. All rights reserved.
////
//
//#include "functions.h"
//#include <fstream>
//#include <iomanip>
//using namespace std;
////#include <gsl/gsl_sf_gamma.h>
////#include <gsl/gsl_sf_psi.h>
////
////gsl_integration_cquad_workspace* workspace;
////
//typedef complex<long double> cpx;
//double clenshaw(long double t);
//cpx ISFs(cpx s, long double* para);
//void NiLT_weeks(cpx (*fun)(cpx, long double*), long double* para);
//fftwStruct iLTStruct;
//
//int main()
//{
//    iLTStruct.fftwIn=fftwl_alloc_complex(M);
//    iLTStruct.fftwOut=fftwl_alloc_complex(M);
//    iLTStruct.integration=fftwl_plan_dft_1d(M, iLTStruct.fftwIn, iLTStruct.fftwOut, FFTW_FORWARD, FFTW_MEASURE);
//    iLTStruct.CoeA.resize(M);
//	long double q=3.0l;
//    long double v0=13.0l;//gsl_vector_get(para, 1);
//    long double lambda=1.0l;
//    long double D=0.4l;//gsl_vector_get(para, 3);\
//
//    long double kv0=q*v0;
//    long double Dq2=D*q*q;
//    long double paraISF[4]={kv0, lambda, Dq2, q};
//    NiLT_weeks(ISFs, paraISF);
//    ofstream debug("debug.txt");
//  
//    //for (int iter = 0; iter<M; ++iter)
//    //{
//    //  debug << setprecision(30) << iLTStruct.CoeA[iter] << endl;
//    //  //Temperary variables used for acceleration.
//    //}
//    for (int iter = 0; iter<num_fit; ++iter)
//    {
//        long double t=(iter+1)*0.01l;
//        double rtd=clenshaw(t);
//        debug << setprecision(30) << rtd << endl;
//        //Temperary variables used for acceleration.
//    }
//
//    fftwl_free(iLTStruct.fftwIn);
//    fftwl_free(iLTStruct.fftwOut);
//    return 0;
//}