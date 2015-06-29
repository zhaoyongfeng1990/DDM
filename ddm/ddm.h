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

#ifdef ISFRunAndTumbleAndDiffusionNoLT
#include "NILT.h"
#endif

#ifdef ISFRTDPNoLT
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
    //For storing the time difference of the imageSeqk
    gsl_matrix* datag;      //g(q,t) matrix.
    gsl_matrix* ldatag;      //g(q,t) matrix.
    
    int qsize;				//Element number of q array.
    vector<double> qabs;	//Absolute value of q array.
    
    double tau[num_fit];
    double s[num_fit];
    
    double inipara[numOfPara];
    gsl_matrix* fittedPara;	//To store the fitting result and error.
    gsl_matrix* fitErr;
    int* status;		//Record the status of fitting.
    gsl_matrix* datafit;
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
    
} dataStruct;

//ISF and its Jacobian. Used in GSL fitting algorithm.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y);
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J);

#endif
