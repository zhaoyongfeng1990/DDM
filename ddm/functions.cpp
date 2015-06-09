//
//  functions.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 14/11/3.
//  Copyright (c) 2014年 ZYF. All rights reserved.
//

#include "functions.h"

//Reading tif files
gsl_matrix* readTiff(const string tifName)
{
    int width=dim;	//Width and height of the image equals to dim (e.g. =512 pixel)
    int height=dim;
    int scale;		//Maximum grayscale, =65535 for 16 bits tif
    TIFFILEHEAD fileHead;	//File head
    TIFFIFD* firstIFD;
    
    ifstream infile(tifName.c_str(),ios_base::binary);
    infile.read((char *)&fileHead.biOrder,2);
    infile.read((char *)&fileHead.version,2);
    infile.read((char *)&fileHead.offsetIFD,4);
    
    firstIFD=new TIFFIFD();
    TIFFIFD* currentIFD=firstIFD;
    infile.seekg(fileHead.offsetIFD,ios::beg);
    uint32_t dataoffset=0;
    int bit=16;
    
    while(1)		//Read the image meta-information from the data structure
    {
        infile.read((char *)&(currentIFD->DEC),2);
        currentIFD->DE=new TIFFDE[currentIFD->DEC];
        for (int c=0; c<currentIFD->DEC; ++c)
        {
            infile.read((char *)&(currentIFD->DE[c].tag),2);
            infile.read((char *)&(currentIFD->DE[c].type),2);
            infile.read((char *)&(currentIFD->DE[c].length),4);
            infile.read((char *)&(currentIFD->DE[c].valueOffset),4);
            if (currentIFD->DE[c].tag==257)
                height=currentIFD->DE[c].valueOffset;
            else if (currentIFD->DE[c].tag==256)
                width=currentIFD->DE[c].valueOffset;
            else if (currentIFD->DE[c].tag==273)
                dataoffset=currentIFD->DE[c].valueOffset;
            else if (currentIFD->DE[c].tag==258)
                bit=currentIFD->DE[c].valueOffset;
        }
        //        if (currentIFD->DE[0].type==1)
        //            scale=256;
        //        else if (currentIFD->DE[0].type==3)
        scale=65536;	//It's often not working to read from meta-information. Set it directly.
        //        else if (currentIFD->DE[0].type==4)
        //            scale=4294967296;
        uint32_t offset;
        infile.read((char *)&(offset),4);
        if (offset==0)
        {
            currentIFD->nextIFD=NULL;
            break;
        }
        else	//IFD of TIFF is a linked list
        {
            infile.seekg(offset,ios::beg);
            currentIFD->nextIFD=new TIFFIFD();
            currentIFD=currentIFD->nextIFD;
        }
    };
    
    gsl_matrix* data=gsl_matrix_calloc(height, width);	//Got the memory for storing image data
    
    //cout << sizeof(TIFFILEHEAD);
    for (int c=height-1; c>-1; --c)	//Reading data. The data of TIFF is stored backward
    {
        infile.seekg(dataoffset/*sizeof(TIFFILEHEAD)*/+(c-height)*width*2,ios::beg);
        for (int c2=0; c2<width; ++c2)
        {
            uint16_t temp=0;
            infile.read((char *)&temp, bit/8);
            gsl_matrix_set(data, c, c2, (double)temp);
        }
    }
    
    //IFD is not needed after we got the data.
    TIFFIFD* next=firstIFD->nextIFD;
    TIFFIFD* current=firstIFD;
    delete current->DE;
    delete current;
    while (next!=NULL)
    {
        current=next;
        next=current->nextIFD;
        delete current->DE;
        delete current;
    };
    
    infile.close();
    return data;
}

gsl_matrix* readSim(const string fileName)	//Read simulation data from ASCII file.
{
    gsl_matrix* data=gsl_matrix_calloc(dim, dim);
    ifstream infile(fileName.c_str());
    for (int iterx=0; iterx<dim; ++iterx)
    {
        for (int itery=0; itery<dim; ++itery)
        {
            double temp;
            infile >> temp;
            gsl_matrix_set(data, iterx, itery, temp);
        }
    }
    return data;
}

///////////////Functions not in used anymore/////////////////////////////////
//The image is shifted by minimizing the correlation between two successive images.
void shiftImage(gsl_matrix** img, const int drow, const int dcol)
{
    gsl_matrix* shifted=gsl_matrix_alloc((*img)->size1, (*img)->size2);
//#pragma omp parallel for
    for (int iterx=0; iterx<(*img)->size1; ++iterx)
    {
        for (int itery=0; itery<(*img)->size2; ++itery)
        {
            gsl_matrix_set(shifted, iterx, itery, gsl_matrix_get((*img), (iterx+drow+dim)%dim, (itery+dcol+dim)%dim));
        }
    }
    gsl_matrix_free((*img));
    (*img)=shifted;
}

//The correlation of shifted two images. The displacement range is [-drow, drow]x[-dcol, dcol]. aveVec is a constant vector. I hope BLAS can speed up the calculation of averaging.
double shiftCorrelation2d(const gsl_matrix* img1, const gsl_matrix* img2, const int drow, const int dcol, const gsl_vector* aveVec)
{
    gsl_vector* temp=gsl_vector_alloc(dim);
    gsl_blas_dgemv(CblasNoTrans, 1, img1, aveVec, 0, temp);
    double meanImg1=0;
    gsl_blas_ddot(temp, aveVec, &meanImg1);
    
    gsl_blas_dgemv(CblasNoTrans, 1, img2, aveVec, 0, temp);
    double meanImg2=0;
    gsl_blas_ddot(temp, aveVec, &meanImg2);
    
    double moment=0;
    double normA=0;
    double normB=0;
//#pragma omp parallel for reduction(+: moment, normA, normB)
    for(int iterrow=0; iterrow<dim; ++iterrow)
        for (int itercol=0; itercol<dim; ++itercol)
        {
            double fluA=gsl_matrix_get(img1, iterrow, itercol)-meanImg1;
            double fluB=gsl_matrix_get(img2, (iterrow+drow+dim)%dim, (itercol+dcol+dim)%dim)-meanImg2;
            moment+=fluA*fluB;
            normA+=fluA*fluA;
            normB+=fluB*fluB;
        }
    
    return moment/sqrt(normA*normB);
}

//The correlation of two images. aveVec is a constant vector. I hope BLAS can speed up the calculation of averaging.
double correlation2d(const gsl_matrix* img1, const gsl_matrix* img2, const gsl_vector* aveVec)
{
    gsl_vector* temp=gsl_vector_alloc(dim);
    gsl_blas_dgemv(CblasNoTrans, 1, img1, aveVec, 0, temp);
    double meanImg1=0;
    gsl_blas_ddot(temp, aveVec, &meanImg1);
    
    gsl_blas_dgemv(CblasNoTrans, 1, img2, aveVec, 0, temp);
    double meanImg2=0;
    gsl_blas_ddot(temp, aveVec, &meanImg2);
    
    double moment=0;
    double normA=0;
    double normB=0;
//#pragma omp parallel for reduction(+: moment, normA, normB)
    for(int iterrow=0; iterrow<dim; ++iterrow)
        for (int itercol=0; itercol<dim; ++itercol)
        {
            double fluA=gsl_matrix_get(img1, iterrow, itercol)-meanImg1;
            double fluB=gsl_matrix_get(img2, iterrow, itercol)-meanImg2;
            moment+=fluA*fluB;
            normA+=fluA*fluA;
            normB+=fluB*fluB;
        }
    
    return moment/sqrt(normA*normB);
}

//Find the position of a particular member in a sorted vector with increasing order, using dichotomy.
int quickFind(const vector<int>& vec, const int value)
{
    int left=0;
    int right=(int)vec.size();
    while(right-left>1)
    {
        int middle=(left+right)/2;
        if (vec[middle]==value)
        {
            return middle;
        }
        else if (vec[middle]>value)
        {
            right=middle;
        }
        else
        {
            left=middle;
        }
    }
    if (vec[left]==value)
        return left;
    else if (vec[right]==value)
        return right;
    else
        return -1;
}

//Find the position of a particular member in a vector.
int find(const vector<int>& vec, const int value)
{
    for (int iter=0; iter<vec.size(); ++iter)
    {
        if (vec[iter]==value)
            return iter;
    }
    return -1;
}
/////////////////////////////////Not in used functions ended////////////////////////

//API for GSL solver.
int fdISFfun(const gsl_vector* para, void* sdata, gsl_vector* y, gsl_matrix* J)
{
    ISFfun(para, sdata, y);
    dISFfun(para, sdata, J);
    
    return  GSL_SUCCESS;
}

//Test fitting result using 0-norm (maximum absolute value) of the relative step size dx/x. tol is the relative tolerance of error, tole is the absolute tolerance of error.
int norm0_rel_test(const gsl_vector * dx, const gsl_vector * x, double tol, double tole)
{
    //Check if some parameters become NAN or INF.
	for (int iter = 0; iter < numOfPara; ++iter)
	{
		if (!isfinite(x->data[iter]))
		{
			return GSL_EOVRFLW;
		}
	}
	for (int iter = 0; iter < numOfPara; ++iter)
	{
		//double relerr = abs(dx->data[iter]/* / x->data[iter]*/);
        if (abs(dx->data[iter]) > tol*abs(x->data[iter])+tole)
		{
			return GSL_CONTINUE;
		}
	}
	return GSL_SUCCESS;
}

//Test fitting result using error estimation from covariance matrix, not reliable. tol is the relative tolerance of error.
int covar_rel_test(const gsl_matrix* J, const gsl_vector* x, double tol)
{
    //Check if some parameters become NAN or INF.
	for (int iter = 0; iter < numOfPara; ++iter)
	{
		if (!isfinite(x->data[iter]))
		{
			return GSL_EOVRFLW;
		}
	}
	gsl_matrix* covar = gsl_matrix_alloc(numOfPara, numOfPara);
	gsl_multifit_covar(J, 0.0, covar);	//Get the covariance matrix.
	double fitErr[numOfPara];
	for (int iterpara = 0; iterpara<numOfPara; ++iterpara)	//err_i=\sqrt{c_{ii}}
	{
		fitErr[iterpara]=sqrt(gsl_matrix_get(covar, iterpara, iterpara));
	}
	gsl_matrix_free(covar);
	for (int iter = 0; iter < numOfPara; ++iter)
	{
		double relerr = abs(fitErr[iter]) / x->data[iter];
		if (relerr>tol)
		{
			return GSL_CONTINUE;
		}
	}
	return GSL_SUCCESS;
}