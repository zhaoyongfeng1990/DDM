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
//#define ISFRUNANDTUMBLE
//#define ISFSWIMMER
//#define ISFSWIMMERSIMPLE
//#define ISFRUNANDTUMBLE_3D
//#define ISFRunAndTumbleAndDiffusion
//#define ISFRunAndTumbleAndDiffusionAndPv
//#define ISFRunAndTumbleAndDiffusionNoLT
//#define ISFRTDPNoLT
#define ISFRTDPNoLT_sigma

#include <cmath>

using namespace std;

const int OMP_NUM_THREADS=1;

const long double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068l;
const double sqrtpi = sqrt(pi); //Constant for convenience
const double s2pi=sqrt(2*pi);

const int dimy = 512;    //Size of the image
const int dimx = 512;
const int dimkx = dimx / 2 + 1;   //Number of wavenumber
const int dimky = dimy / 2 + 1;   //Number of wavenumber
const int numOfSeq = 4500;  //Number of total time points in experiment
const int numOfDiff = 4000; //Number of tau
const int numOfk = dimy*dimkx;    //Number of data points after FFTpoints used in fitting.

const double dx = 6.5 / 4.0; // 0.65;   //Pixel size
const double dqx = 2 * pi / dx / dimx;    //q step after FFT
const double dqy = 2 * pi / dx / dimy;    //q step after FFT
const double qmax=2 * pi / dx /sqrt(2);   //Maximum possible value of q
const double qstep=0.01;    //The width of cirque when averaging the direction of q.

const double dt = 1.0 / 100;    //Time step
const double ds=0.01;
const double smin=0.01;


const int maxIter=200;    //Maximum iteration number in fitting

const double alphaGuess=0.8;
const double DGuess=0.4;//0.2;
const double vbarGuess=13;
const double lambdaGuess=1;
const double ZGuess=3;
const double sigmaGuess=vbarGuess/8;

//const int winDim = 5;   //Size of searching window. Used in aligning images.

//Number of model parameters
#ifdef ISFSWIMMER
const int numOfPara = 6;
//#define NoJacobian
#endif

#ifdef ISFSWIMMERSIMPLE
const int numOfPara = 5;
#endif

#ifdef ISFRUNANDTUMBLE
const int numOfPara = 4;
const double precision = 1e-15; //Used in 2D R&T model. The precision of numerical evaluation.
#endif

#ifdef ISFRUNANDTUMBLE_3D
const int numOfPara=4;
#define NeedLaplaceTrans
#endif

#ifdef ISFRunAndTumbleAndDiffusion
const int numOfPara=6;
#define NeedLaplaceTrans
#endif

#ifdef ISFRunAndTumbleAndDiffusionNoLT
const int numOfPara=6;
//#define NoJacobian
#define NeedNumericalInverseLaplaceTransformation
#endif

#ifdef ISFRunAndTumbleAndDiffusionAndPv
const int numOfPara=7;
#define NeedLaplaceTrans
#endif

#ifdef ISFRTDPNoLT
const int numOfPara=7;
//#define NoJacobian
#define NeedNumericalInverseLaplaceTransformation
#define IfComplexIntegration
#endif

#ifdef ISFRTDPNoLT_sigma
const int numOfPara=7;
//#define NoJacobian
#define NeedNumericalInverseLaplaceTransformation
#define IfComplexIntegration
#define MultiQFit
#endif

#endif /* defined(__ddm__parameters__) */
