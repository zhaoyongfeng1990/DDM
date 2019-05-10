//
//  ddm.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <fstream>
#include <omp.h>

//Initialization. Calculation of some constant.
ddm_Multisets::ddm_Multisets() : qabs_low(), qabs_high(), tau_low(), tau_high()
{
	datag_low = nullptr;
	datag_high = nullptr;
    fittedPara=nullptr;
    fitErr=nullptr;
    status=nullptr;
    qsize_low=0;
	qsize_high = 0;
    
    ifstream paraFile("parameters.txt");
    
    paraFile >> OMP_NUM_THREADS;
	paraFile >> num_qCurve_low;
	qIncreList_low = new int[num_qCurve_low];
	qIncreList_low[0] = 0;
	for (int i = 1; i<num_qCurve_low; ++i)
	{
		paraFile >> qIncreList_low[i];
	}
	paraFile >> num_qCurve_high;
	qIncreList_high = new int[num_qCurve_high];
	qIncreList_high[0] = 0;
	for (int i = 1; i<num_qCurve_high; ++i)
	{
		paraFile >> qIncreList_high[i];
	}

	paraFile >> timeWindow_low;
	paraFile >> timeWindow_high;
    
    paraFile >> maxIter;
    
    paraFile >> alphaGuess;
    paraFile >> DGuess;
    paraFile >> vbarGuess;
    paraFile >> lambdaGuess;
    paraFile >> sigmaGuess;
    paraFile >> TTGuess;
    
    paraFile.close();
    
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
ddm_Multisets::~ddm_Multisets()
{
	if (datag_low != nullptr)
		gsl_matrix_free(datag_low);
	if (datag_high != nullptr)
		gsl_matrix_free(datag_high);
    if (fittedPara!=nullptr)
        gsl_matrix_free(fittedPara);
    if (fitErr!=nullptr)
        gsl_matrix_free(fitErr);
    if (status!=nullptr)
        delete [] status;
	if (qIncreList_low != nullptr)
	{
		delete[] qIncreList_low;
	}
	if (qIncreList_high != nullptr)
	{
		delete[] qIncreList_high;
	}
}