////
////  debugmain.cpp
////  ddm
////
////  Created by Zhao Yongfeng on 15/6/5.
////  Copyright (c) 2015年 ZYF. All rights reserved.
////
//
//#include "functions.h"
//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_sf_psi.h>
//
//gsl_integration_cquad_workspace* workspace;
//
//int main()
//{
//#ifdef ISFRunAndTumbleAndDiffusionAndPv
//    workspace=gsl_integration_cquad_workspace_alloc(100000);
//#endif
//    double IntPara[14];
//    double vbar=10;
//    double Z=500;
//    double lambda=1;
//    double D=0.4;
//    double q=0.01;
//    
//    IntPara[0]=vbar;
//    IntPara[1]=Z;
//    IntPara[2]=lambda;
//    
//    double dq2=D*q*q;
//    IntPara[3]=q;
//    double Z1=Z+1;
//    IntPara[5]=Z1/vbar;
//    
//    IntPara[6]=Z1*log(IntPara[5])-gsl_sf_lngamma(Z1);
//    IntPara[7]=Z1/vbar/vbar;
//    
//    IntPara[8]=gsl_sf_psi(Z1);
//    
//    IntPara[9]=pow(q,4);
//    gsl_function Fun;
//    Fun.function=&integrandFun;
//    Fun.params=IntPara;
//    
//    double tau=0.4;
//    IntPara[4]=tau+lambda+dq2;
//    IntPara[10]=IntPara[4]*IntPara[4];
//    IntPara[13]=0;
//    
//    double pts[2]={0,1};
//    double yi=0;
//    double err=0;
//    size_t nevals;
//    gsl_integration_cquad(&Fun, 0, 1, 0, 1e-7, workspace, &yi, &err, &nevals);
//    cout << yi << "+-" << err << endl;
//    
//#ifdef ISFRunAndTumbleAndDiffusionAndPv
//        gsl_integration_cquad_workspace_free(workspace);
//#endif
//    return 0;
//}