//
//  shiftImage.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

///////////////////////////////////////////////////////
//Functions for shifting images, not in use any more!//
///////////////////////////////////////////////////////

#include "ddm.h"
#include <omp.h>

//void shiftImage(gsl_matrix** img, const int drow, const int dcol);      //Used for aligning viberating images.
//double shiftCorrelation2d(const gsl_matrix* img1, const gsl_matrix* img2, const int drow, const int dcol, const gsl_vector* aveVec);
////Calculating shifted correlation for two images, drow and dcol indicate the displacement between two images.
//double correlation2d(const gsl_matrix* img1, const gsl_matrix* img2, const gsl_vector* aveVec);
//Calculating correlation for two images.


void ddm::shiftImages()
{
//    for (int iter = 0; iter < numOfDiff; ++iter)
//    {
//    	double max = 0;
//    	int drow = 0;
//    	int dcol = 0;
//    	for (int iterrow = -winDim; iterrow <= winDim; ++iterrow)
//    	{
//    		for (int itercol = -winDim; itercol <= winDim; ++itercol)
//    		{
//    			double corr = shiftCorrelation2d(imageSeq[iter], imageSeq[iter + 1], iterrow, itercol, aveVec);
//    			if (corr > max)
//    			{
//    				max = corr;
//    				drow = iterrow;
//    				dcol = itercol;
//    			}
//    		}
//    	}
//    	shiftImage(&imageSeq[iter + 1], drow, dcol);
//    }
}


////The image is shifted by minimizing the correlation between two successive images.
//void shiftImage(gsl_matrix** img, const int drow, const int dcol)
//{
//    gsl_matrix* shifted=gsl_matrix_alloc((*img)->size1, (*img)->size2);
//#pragma omp parallel for
//    for (int iterx=0; iterx<(*img)->size1; ++iterx)
//    {
//        for (int itery=0; itery<(*img)->size2; ++itery)
//        {
//            gsl_matrix_set(shifted, iterx, itery, gsl_matrix_get((*img), (iterx+drow+dimx)%dimx, (itery+dcol+dimy)%dimy));
//        }
//    }
//    gsl_matrix_free((*img));
//    (*img)=shifted;
//}
//
////The correlation of shifted two images. The displacement range is [-drow, drow]x[-dcol, dcol]. aveVec is a constant vector. I hope BLAS can speed up the calculation of averaging.
//double shiftCorrelation2d(const gsl_matrix* img1, const gsl_matrix* img2, const int drow, const int dcol, const gsl_vector* aveVec)
//{
//    gsl_vector* temp=gsl_vector_alloc(dim);
//    gsl_blas_dgemv(CblasNoTrans, 1, img1, aveVec, 0, temp);
//    double meanImg1=0;
//    gsl_blas_ddot(temp, aveVec, &meanImg1);
//    
//    gsl_blas_dgemv(CblasNoTrans, 1, img2, aveVec, 0, temp);
//    double meanImg2=0;
//    gsl_blas_ddot(temp, aveVec, &meanImg2);
//    
//    double moment=0;
//    double normA=0;
//    double normB=0;
//#pragma omp parallel for reduction(+: moment, normA, normB)
//    for(int iterrow=0; iterrow<dimy; ++iterrow)
//        for (int itercol=0; itercol<dimx; ++itercol)
//        {
//            double fluA=gsl_matrix_get(img1, iterrow, itercol)-meanImg1;
//            double fluB=gsl_matrix_get(img2, (iterrow+drow+dimy)%dimy, (itercol+dcol+dimx)%dimx)-meanImg2;
//            moment+=fluA*fluB;
//            normA+=fluA*fluA;
//            normB+=fluB*fluB;
//        }
//    
//    return moment/sqrt(normA*normB);
//}
//
////The correlation of two images. aveVec is a constant vector. I hope BLAS can speed up the calculation of averaging.
//double correlation2d(const gsl_matrix* img1, const gsl_matrix* img2, const gsl_vector* aveVec)
//{
//    gsl_vector* temp=gsl_vector_alloc(dim);
//    gsl_blas_dgemv(CblasNoTrans, 1, img1, aveVec, 0, temp);
//    double meanImg1=0;
//    gsl_blas_ddot(temp, aveVec, &meanImg1);
//    
//    gsl_blas_dgemv(CblasNoTrans, 1, img2, aveVec, 0, temp);
//    double meanImg2=0;
//    gsl_blas_ddot(temp, aveVec, &meanImg2);
//    
//    double moment=0;
//    double normA=0;
//    double normB=0;
//#pragma omp parallel for reduction(+: moment, normA, normB)
//    for(int iterrow=0; iterrow<dim; ++iterrow)
//        for (int itercol=0; itercol<dim; ++itercol)
//        {
//            double fluA=gsl_matrix_get(img1, iterrow, itercol)-meanImg1;
//            double fluB=gsl_matrix_get(img2, iterrow, itercol)-meanImg2;
//            moment+=fluA*fluB;
//            normA+=fluA*fluA;
//            normB+=fluB*fluB;
//        }
//    
//    return moment/sqrt(normA*normB);
//}
