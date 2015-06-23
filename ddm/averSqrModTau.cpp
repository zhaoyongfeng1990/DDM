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
    int progress = 0;       //Indicator of the progess.
#pragma omp parallel for
    for (int iterdiff = 1; iterdiff <= numOfDiff; ++iterdiff)
    {
        gsl_matrix* temp = gsl_matrix_alloc(dimy, dimkx);
        gsl_matrix_set_zero(temp);
        for (int itert = 0; itert < numOfSeq - numOfDiff; ++itert)
        {
            for (int itermem = 0; itermem < numOfk; ++itermem)
            {
                double* later = imageSeqk[itert + iterdiff]->data;		//I(q, t+\tau)
                double* before = imageSeqk[itert]->data;				//I(q, t)
                double real = later[itermem * 2] - before[itermem * 2];	//Difference of Re
                double image = later[itermem * 2 + 1] - before[itermem * 2 + 1];	//Difference of Im
                temp->data[itermem] += real*real + image*image;			//Noticed |a|=|a*|
            }
        }
        progress += 1; // numOfSeq - iterdiff;
        gsl_matrix_scale(temp, 1.0 / (numOfSeq - numOfDiff));		//Average on t
        imagekDiff[iterdiff - 1] = temp;
        cout << "Calculating average of square module for different tau... " << 100.0*progress / numOfDiff << "% finished." << endl;
    }
}