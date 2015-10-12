//
//  parameters.h
//  ddm
//
//  Created by Zhao Yongfeng on 14/11/3.
//  Copyright (c) 2014年 ZYF. All rights reserved.
//

#ifndef __ddm__parameters__
#define __ddm__parameters__

//If the platform is Windows, uncomment this line.
//#define WINDOWS
#ifdef WINDOWS
#pragma comment(lib,"libgsl.lib")
#pragma comment(lib, "libgslcblas.lib")
#pragma comment(lib, "libfftw3-3.lib")
#include <cstdint>
#endif
///////////////////

//Switches of different models
//#include "ISFtype.h"
//#define ISFRDP
//#define ISFRD
//#define ISFRTD
//#define ISFRTDP
#define ISFRTDPTT

#include <cmath>

using namespace std;


const long double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068l;

//If you found some problem with the Jaccobian calculation, you can try to turn on this switch.
//#define NoJacobian

//Number of model parameters
#ifdef ISFRDP
const int numOfPara = 4;
#endif

#ifdef ISFRD
const int numOfPara = 3;
#endif

#ifdef ISFRTD
const int numOfPara=4;
#define NeedNumericalInverseLaplaceTransformation
#endif

#ifdef ISFRTDP
const int numOfPara=5;
#define NeedNumericalInverseLaplaceTransformation
#define IfComplexIntegration
#endif

#ifdef ISFRTDPfix
const int numOfPara=1;
#define NeedNumericalInverseLaplaceTransformation
#define IfComplexIntegration
#endif

#ifdef ISFRTDPTT
const int numOfPara=6;
#define NeedNumericalInverseLaplaceTransformation
#define IfComplexIntegration
#endif

#endif /* defined(__ddm__parameters__) */
