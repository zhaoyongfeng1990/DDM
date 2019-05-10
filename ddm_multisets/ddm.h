//
//  ddm.h
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef ddm_ddm_h
#define ddm_ddm_h

#include <gsl/gsl_matrix.h>
#include <vector>
#include <fftw3.h>	//FFTW
#include <string>
#include "parameters.h"

#ifdef NeedNumericalInverseLaplaceTransformation
#include "NILT.h"
#endif

class ddm_Multisets
{
public:
    //Methods
    ddm_Multisets();
    ~ddm_Multisets();
    
    //Reading g(q,t), q, and tau data from corresponding files.
    void recover();
    
    //Fitting.
    void fitting();
    //Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.
    int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole);
    
    //Print fitting result
    void printFit();
    //Print arbitrary gsl_matrix m for debug
    void printdebugM(gsl_matrix* m, const string filename);
    
    /////////////////////////////////////////////////////
    //Variables
	vector<double> qabs_low;	//Absolute value of q
	vector<double> qabs_high;	//Absolute value of q
	vector<double> tau_low;     //List of time point
	vector<double> tau_high;     //List of time point
    //For storing the time difference of the imageSeqk
	gsl_matrix* datag_low;      //g(q,t) matrix
	gsl_matrix* datag_high;      //g(q,t) matrix
    gsl_matrix* fittedPara;	//To store the fitting result
    gsl_matrix* fitErr;     //To store the fitting error
	double inipara[numOfPara];  //Initial parameters
    int* status;		//Record the status of fitting
	int* qIncreList_low;    //The increment of q for multiple q curves
	int* qIncreList_high;    //The increment of q for multiple q curves
    
    int qsize_low;	 //Element number of q array.
	int qsize_high;
    
    int OMP_NUM_THREADS;    //Number of threads for openMP
	int num_fit_low;   //Number of data points used in fitting.
	int num_fit_high;   //Number of data points used in fitting.
	int num_qCurve_low; //Number of q curves
	int num_qCurve_high; //Number of q curves
	double timeWindow_low;  //Time window that is used in fitting
	double timeWindow_high;  //Time window that is used in fitting
	double maxIter;    //Maximum iteration number in fitting

					   //The initial values to be passed from outside
	double alphaGuess;
	double DGuess;
	double vbarGuess;
	double lambdaGuess;
	double sigmaGuess;
	double TTGuess;
};

//Data stuct used in GSL fitting algorithm.
typedef struct
{
	double* data;
	double* q;
	double* tau;
	int num_fit;
	int num_qCurve_low;
	int num_qCurve_high;
    
#ifdef ISFRTD
    NILT* ISFILT;
    NILT* dvISFILT;
    NILT* dDISFILT;
    NILT* dlambdaISFILT;
#endif
    
#ifdef ISFRTDP
    NILT* ISFILT;
    NILT* dvbarISFILT;
    NILT* dsigmaISFILT;
    NILT* dDISFILT;
    NILT* dlambdaISFILT;
#endif
    
#ifdef ISFRTDPTT
    NILT* ISFILT;
    NILT* dvbarISFILT;
    NILT* dsigmaISFILT;
    NILT* dDISFILT;
    NILT* dlambdaISFILT;
    NILT* dTTISFILT;
#endif
    
#ifdef ISFRTDPfix
    long double alpha;
    long double D;
    long double vbar;
    long double sigma;
    
    long double vbsigma2;
    long double logfactor;
    long double vb2sigma2;
    long double cpsiz1;
    long double vb2sigma3;
    
    NILT* ISFILT;
    NILT* dlambdaISFILT;
#endif
    
#ifdef ISFRTDPTTfix
    long double alpha;
    long double D;
    long double vbar;
    long double sigma;
    
    long double vbsigma2;
    long double logfactor;
    long double vb2sigma2;
    long double cpsiz1;
    long double vb2sigma3;
    
    NILT* ISFILT;
    NILT* dlambdaISFILT;
    NILT* dTTISFILT;
#endif
    
} dataStruct;

//ISF and its Jacobian. Used in GSL fitting algorithm.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y);
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J);

#endif
