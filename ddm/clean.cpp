//
//  clean.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"

//Free the memory used by imageSeqk
void ddm::cleanSeqk()
{
    size_t size=imageSeqk.size();
#pragma omp parallel for
    for (int iter = 0; iter < size; ++iter)		//Free the memory
        gsl_matrix_complex_free(imageSeqk[iter]);
    imageSeqk.clear();
}

//Free the memory used by imagekDiff
void ddm::cleankDiff()
{
    size_t size=imagekDiff.size();
#pragma omp parallel for
    for (int iter = 0; iter < size; ++iter)
        gsl_matrix_free(imagekDiff[iter]);
    imagekDiff.clear();
}