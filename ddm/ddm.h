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
    
    void recover();
    
    void readAndFFT(const string filePrefix);
    void shiftImages();
    
    void averSqrModTau();
    
    void aveQSort();
    void aveQConcentric();
    void aveQBilinear();
    
    void LaplaceTrans();
    
    void fitting();
    void fitting_estRange();
    
    void printG();
    void printGs();
    void printFit();
    
    void cleanSeqk();
    void cleankDiff();
    
    /////////////////////////////////////////////////////
    //Variables
    
    //vector<gsl_matrix*> imageSeq(numOfSeq);
    ////Sequence of images
    vector<gsl_matrix_complex*> imageSeqk;
    //Sequence for storing image after FFT.
    vector<gsl_matrix*> imagekDiff;
    vector<double> qabs;	//Absolute value of q array.
    //For storing the time difference of the imageSeqk
    gsl_matrix* datag;      //g(q,t) matrix.
    gsl_matrix* ldatag;      //g(q,t) matrix.
    gsl_matrix* datafit;
    gsl_matrix* fittedPara;	//To store the fitting result and error.
    gsl_matrix* fitErr;
    int* status;		//Record the status of fitting.
    double tau[num_fit];
    double s[num_fit];
    double inipara[numOfPara];
    
    int qsize;				//Element number of q array.
    
    int iniTime;
    int finalTime;
    int num_fit;
    //gsl_vector* aveVec;
};

//Data stuct used in GSL fitting algorithm.
typedef struct
{
    double* data;
    double* tau;
    double q;
    
#ifdef ISFRunAndTumbleAndDiffusionNoLT
    NILT* ISFILT;
    NILT* dvISFILT;
    NILT* dDISFILT;
    NILT* dlambdaISFILT;
#endif
    
#ifdef ISFRTDPNoLT
    NILT* ISFILT;
    NILT* dvbarISFILT;
    NILT* dZISFILT;
    NILT* dDISFILT;
    NILT* dlambdaISFILT;
#endif
    
#ifdef ISFRTDPNoLT_sigma
    NILT* ISFILT;
    NILT* dvbarISFILT;
    NILT* dsigmaISFILT;
    NILT* dDISFILT;
    NILT* dlambdaISFILT;
#endif
    
} dataStruct;

//ISF and its Jacobian. Used in GSL fitting algorithm.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y);
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J);

#endif
