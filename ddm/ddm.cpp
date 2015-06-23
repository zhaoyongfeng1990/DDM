//
//  ddm.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <fftw3.h>	//FFTW
#include "ddm.h"

ddm::ddm() : imageSeqk(numOfSeq), imagekDiff(numOfDiff), qabs()
{
    datag=nullptr;
    ldatag=nullptr;
    fittedPara=nullptr;
    fitErr=nullptr;
    status=nullptr;
    datafit=nullptr;
    qsize=0;
    
    for (int itertau=0; itertau<num_fit; ++itertau)
    {
#ifdef NeedLaplaceTrans
        //s[itertau]=exp(smin+ds*itertau);
        s[itertau]=smin+ds*itertau;
#else
        tau[itertau]=(itertau+1)*dt;
#endif
    }
    
#ifdef ISFSWIMMER
    inipara[0]=alphaGuess;
    inipara[1]=DGuess;
    inipara[2]=vbarGuess;
    inipara[3]=ZGuess;
#endif
    
#ifdef ISFSWIMMERSIMPLE
    inipara[0]=alphaGuess;
    inipara[1]=DGuess;
    inipara[2]=vbarGuess;
#endif
    
#ifdef ISFRUNANDTUMBLE
    inipara[0]=lambdaGuess;
    inipara[1]=vbarGuess;
#endif
    
#ifdef ISFRUNANDTUMBLE_3D
    inipara[0]=vbarGuess;
    inipara[1]=lambdaGuess;
#endif
    
#ifdef ISFRunAndTumbleAndDiffusion
    inipara[0]=alphaGuess;
    inipara[1]=vbarGuess;
    inipara[2]=lambdaGuess;
    inipara[3]=DGuess;
#endif
    
#ifdef ISFRunAndTumbleAndDiffusionAndPv
    inipara[0]=alphaGuess;
    inipara[1]=vbarGuess;
    inipara[2]=ZGuess;
    inipara[3]=lambdaGuess;
    inipara[4]=DGuess;
#endif

#ifdef ISFRunAndTumbleAndDiffusionNoLT
    inipara[0]=alphaGuess;
    inipara[1]=vbarGuess;
    inipara[2]=lambdaGuess;
    inipara[3]=DGuess;
#endif
    
    //aveVec = gsl_vector_alloc(dim);
    //gsl_vector_set_all(aveVec, 1.0 / dim);			//Vecter used in calculating average. I hope BLAS can help speed up the average.
    
}

ddm::~ddm()
{
    if (datag!=nullptr)
        gsl_matrix_free(datag);
    if (fittedPara!=nullptr)
        gsl_matrix_free(fittedPara);
    if (fitErr!=nullptr)
        gsl_matrix_free(fitErr);
    if (datafit!=nullptr)
        gsl_matrix_free(datafit);
    if (status!=nullptr)
        delete[] status;
    
    if (imagekDiff.size()!=0)
    {
        cleankDiff();
    }
    if (imageSeqk.size()!=0)
    {
        cleanSeqk();
    }
    
#ifdef NeedLaplaceTrans
    if (ldatag!=nullptr)
        gsl_matrix_free(ldatag);
#endif
}