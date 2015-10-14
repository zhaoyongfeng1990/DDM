//
//  ddm.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <fftw3.h>	//FFTW
#include "ddm.h"
#include <fstream>
#include <omp.h>

//Initialization. Calculation of some constant.
ddm::ddm() : imageSeqk(), imagekDiff(), qabs(), tau()
{
    
    datag=nullptr;
    fittedPara=nullptr;
    fitErr=nullptr;
    status=nullptr;
    qsize=0;
    
    ifstream paraFile("parameters.txt");
    
    paraFile >> OMP_NUM_THREADS;
    paraFile >> num_qCurve;
    qIncreList=new int[num_qCurve];
    qIncreList[0]=0;
    for (int i=1; i<num_qCurve; ++i)
    {
        paraFile >> qIncreList[i];
    }
    
    paraFile >> dimx;
    paraFile >> dimy;
    
    dimkx = dimx / 2 + 1;
    dimky = dimy / 2 + 1;
    numOfk = dimy*dimkx;
    
    paraFile >> numOfSeq;
    paraFile >> numOfDiff;
    paraFile >> dx;
    
    dqx = 2 * pi / dx / dimx;    //q step after FFT
    dqy = 2 * pi / dx / dimy;    //q step after FFT
    qmax=2 * pi / dx /sqrt(2);   //Maximum possible value of q
    
    paraFile >> qmin;
    paraFile >> qstep;
    paraFile >> dt;
    paraFile >> timeWindow;
    
    paraFile >> maxIter;
    
    paraFile >> alphaGuess;
    paraFile >> DGuess;
    paraFile >> vbarGuess;
    paraFile >> lambdaGuess;
    paraFile >> ZGuess;
    paraFile >> sigmaGuess;
    paraFile >> TTGuess;
    
    imageSeqk.reserve(numOfSeq);
    imagekDiff.reserve(numOfDiff);
    
    paraFile.close();
#ifdef ISFRDP
    inipara[0]=alphaGuess;
    inipara[1]=DGuess;
    inipara[2]=vbarGuess;
    inipara[3]=ZGuess;
#endif
    
#ifdef ISFRD
    inipara[0]=alphaGuess;
    inipara[1]=DGuess;
    inipara[2]=vbarGuess;
#endif

#ifdef ISFRTD
    inipara[0]=alphaGuess;
    inipara[1]=vbarGuess;
    inipara[2]=lambdaGuess;
    inipara[3]=DGuess;
#endif
    
#ifdef ISFRTDP
    inipara[0]=alphaGuess;
    inipara[1]=vbarGuess;
    inipara[2]=sigmaGuess;
    inipara[3]=lambdaGuess;
    inipara[4]=DGuess;
#endif
    
#ifdef ISFRTDPfix
    inipara[0]=lambdaGuess;
#endif
    
#ifdef ISFRTDPTT
    inipara[0]=alphaGuess;
    inipara[1]=vbarGuess;
    inipara[2]=sigmaGuess;
    inipara[3]=lambdaGuess;
    inipara[4]=DGuess;
    inipara[5]=TTGuess;
#endif
    
#ifdef ISFRTDPTTfix
    inipara[0]=lambdaGuess;
    inipara[1]=TTGuess;
#endif
    omp_set_num_threads(OMP_NUM_THREADS);     //Set the number of threads
}

//Destructor, cleaning up.
ddm::~ddm()
{
    if (datag!=nullptr)
        gsl_matrix_free(datag);
    if (fittedPara!=nullptr)
        gsl_matrix_free(fittedPara);
    if (fitErr!=nullptr)
        gsl_matrix_free(fitErr);
    if (status!=nullptr)
        delete [] status;
    if (qIncreList!=nullptr)
    {
        delete [] qIncreList;
    }
    if (imagekDiff.size()!=0)
    {
        cleankDiff();
    }
    if (imageSeqk.size()!=0)
    {
        cleanSeqk();
    }
}