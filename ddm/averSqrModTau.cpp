//
//  averSqrModTau.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>
#include <omp.h>

void ddm::averSqrModTau()
{
    //Anti-aliasing
    int cnumOfDiff=numOfDiff;
    int numAve=numOfSeq - cnumOfDiff;
    int cnumOfk=numOfk;
    int cdimkx=dimkx;
    int cdimy=dimy;
    
    //Generate the list of tau, which of the form like n*2^m while n is an odd number given by the initial value of 'candidates'.
    vector<int> diff;   //To store the index of tau value.
    diff.reserve(cnumOfDiff);
    tau.reserve(cnumOfDiff);
    diff.push_back(1);
    tau.push_back(dt);
    int candidates[16]={2, 3, 5, 7, 9,11,13,15,17,19,21,23,25,27,29,31};
    //We always put the first element in candidates to diff, then double this element and do sorting again. As we only have a list of 16 element which is almost sorted, direct insert sorting is good enough.
    while (candidates[0]<cnumOfDiff)
    {
        diff.push_back(candidates[0]);
        tau.push_back(dt*candidates[0]);
        candidates[0]*=2;
        int insertPos=15;
        for (int iteri=1; iteri<16; ++iteri)
        {
            if (candidates[0]<candidates[iteri])
            {
                insertPos=iteri-1;
                break;
            }
        }
        if (insertPos>0)
        {
            int temp=candidates[0];
            for (int iteri=0; iteri<insertPos; ++iteri)
            {
                candidates[iteri]=candidates[iteri+1];
            }
            candidates[insertPos]=temp;
        }
    }
    //We also include the last time point
    diff.push_back(cnumOfDiff);
    tau.push_back(cnumOfDiff*dt);
    num_fit=(int)diff.size();
    
    int cnum_fit=num_fit;
    
    int progress = 0;       //Indicator of the progess.
#pragma omp parallel for
    for (int iterdiff = 0; iterdiff < cnum_fit; ++iterdiff)
    {
        int cdiff=diff[iterdiff];
        gsl_matrix* temp = gsl_matrix_alloc(cdimy, cdimkx);
        gsl_matrix_set_zero(temp);
        for (int itert = 0; itert < numAve; ++itert)
        {
            double* before = imageSeqk[itert]->data;				//I(q, t)
            double* later = imageSeqk[itert + cdiff]->data;		//I(q, t+\tau)
            for (int itermem = 0; itermem < cnumOfk; ++itermem)
            {
                double beforeRe=before[itermem*2];
                double beforeIm=before[itermem*2+1];
                double laterRe=later[itermem*2];
                double laterIm=later[itermem*2+1];
                double real = laterRe - beforeRe;	//Difference of Re
                double image = laterIm - beforeIm;	//Difference of Im
                temp->data[itermem] += real*real + image*image;			//Noticed |a|=|a*|
            }
        }
        progress += 1; // numOfSeq - iterdiff;
        gsl_matrix_scale(temp, 1.0 / numAve);		//Average on t
        imagekDiff[iterdiff] = temp;
        cout << "Calculating average of square module for different tau... " << 100.0*progress / cnum_fit << "% finished." << '\n';
    }
}