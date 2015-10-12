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

class ddm
{
public:
    //Methods
    ddm();
    ~ddm();
    
    //Reading g(q,t), q, and tau data from corresponding files.
    void recover();
    //Reading data from images or simulation result (binary or text), and do FFT on it.
    void readAndFFT(const string filePrefix);
    
    //Calculating <|\Delta I(q, t)|^2>_t
    void averSqrModTau();
    
    //Average the directions of q, by integrating bilinear interpolation function.
    void aveQBilinear();
    
    //Not finished yet!!
    void aveQBicubic();
    
    //Fitting.
    void fitting();
    //Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.
    int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole);
    
    //Print unfitted data
    void printG();
    //Print fitting result
    void printFit();
    //Print arbitrary gsl_matrix m for debug
    void printdebugM(gsl_matrix* m, const string filename);
    
    //Free the memory used by imageSeqk
    void cleanSeqk();
    //Free the memory used by imagekDiff
    void cleankDiff();
    
    /////////////////////////////////////////////////////
    //Variables
    ////Sequence of images
    vector<gsl_matrix_complex*> imageSeqk;
    //Sequence for storing image after FFT
    vector<gsl_matrix*> imagekDiff;
    vector<double> qabs;	//Absolute value of q
    vector<double> tau;     //List of time point
    //For storing the time difference of the imageSeqk
    gsl_matrix* datag;      //g(q,t) matrix
    gsl_matrix* fittedPara;	//To store the fitting result
    gsl_matrix* fitErr;     //To store the fitting error
    int* status;		//Record the status of fitting
    int* qIncreList;    //The increment of q for multiple q curves
    double inipara[numOfPara];  //Initial parameters
    double dx;     //Pixel size
    double dqx;    //q step after FFT
    double dqy;    //q step after FFT
    double qmin;    //Minimal q value
    double qmax;   //Maximum possible value of q
    double qstep;    //The width of cirque when averaging the direction of q.
    double dt;    //Time step
    double maxIter;    //Maximum iteration number in fitting
    
    //The initial values to be passed from outside
    double alphaGuess;
    double DGuess;
    double vbarGuess;
    double lambdaGuess;
    double ZGuess;
    double sigmaGuess;
    
    int qsize;	 //Element number of q array.
    
    int OMP_NUM_THREADS;    //Number of threads for openMP
    int dimy;    //Size of the image
    int dimx;    //Size of the image
    int dimkx;   //Number of wavenumber
    int dimky;   //Number of wavenumber
    int numOfSeq;  //Number of total time points in experiment
    int numOfDiff; //Size of tau array
    int numOfk;    //Number of points in q lattice
    int num_fit;   //Number of data points used in fitting.
    int num_qCurve; //Number of q curves
    double timeWindow;  //Time window that is used in fitting
};

//Data stuct used in GSL fitting algorithm.
typedef struct
{
    double* data;
    double* tau;
    double* q;
    int num_fit;
    int num_qCurve;
    
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
    
} dataStruct;

//ISF and its Jacobian. Used in GSL fitting algorithm.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y);
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J);

#endif
