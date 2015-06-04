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
#pragma omp parallel for
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
#pragma omp parallel for reduction(+: moment, normA, normB)
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
#pragma omp parallel for reduction(+: moment, normA, normB)
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

#ifdef ISFSWIMMER
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters.
    double alpha=gsl_vector_get(para, 0);
    double D=gsl_vector_get(para, 1);
    double v=gsl_vector_get(para, 2);
    double Z=gsl_vector_get(para, 3);
    double A=gsl_vector_get(para, 4);
    double B=gsl_vector_get(para, 5);
    
	for (int iter = 0; iter<num_fit; ++iter)
    {
        double Gamma=q*tau[iter]*v/(Z+1);
        double yi=log(A*(1-exp(-D*q*q*tau[iter])*(1-alpha+alpha/Z/Gamma*sin(Z*atan(Gamma))/pow(1+Gamma*Gamma, Z/2.0)))+B);

		double result = yi - dataAry[iter];

		//Punishment terms, to make constrains in parameter space. 
		if (alpha < 0)
		{
		    result += 1e5 * alpha *alpha;
		}
		
		if (alpha>1)
		{
			result += 1e5*(alpha - 1)*(alpha - 1);
		}

		if (D<0)
		{
			result += 1e5*D*D;
		}

		if (v<0)
		{
			result += 1e5*v*v;
		}

        if (Z<0)
        {
            result += 1e5*Z*Z;
        }
        if (A<0)
        {
            result += 1e5*A*A;
        }

        gsl_vector_set(y, iter, result);
    }

    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters
    double alpha = gsl_vector_get(para, 0);
    double D = gsl_vector_get(para, 1);
    double v = gsl_vector_get(para, 2);
    double Z = gsl_vector_get(para, 3);
	double A = gsl_vector_get(para, 4);
	double B = gsl_vector_get(para, 5);
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double Gamma=q*tau[iter]*v/(Z+1);
        double difexp=exp(-D*q*q*tau[iter]);
        double zatang=Z*atan(Gamma);
        double sinzg=sin(zatang);
        double coszg=cos(zatang);
        double powg=pow(1+Gamma*Gamma, Z/2.0);
        double qqttvv=q*q*tau[iter]*tau[iter]*v*v;
		double zzqqttvv = (1 + Z)*(1 + Z) + qqttvv; 
		double yi = A*(1 - difexp*(1 - alpha + alpha / Z / Gamma*sinzg / powg)) + B;
        //(log y)'=y'/y
        
        gsl_matrix_set(J, iter, 0, (A*difexp*(1-sinzg/Z/Gamma/powg))/yi );
        gsl_matrix_set(J, iter, 1, (A*q*q*tau[iter]*(1-alpha+alpha*sinzg/powg/Z/Gamma)*difexp)/yi );
        
        gsl_matrix_set(J, iter, 2, ((alpha*(1+Z)*A*difexp*((1+Z+qqttvv)*sinzg-v*Z*q*tau[iter]*coszg) )/powg/(v*Z*Gamma*zzqqttvv) )/yi );
        
        gsl_matrix_set(J, iter, 3, ((A*alpha*difexp*( (v*Z*Z*q*tau[iter]-zzqqttvv*zatang)*2*coszg+(2*(1+Z-(Z-1)*qqttvv)+Z*zzqqttvv*log(1+Gamma*Gamma))*sinzg))/powg/Gamma/2/Z/Z/zzqqttvv)/yi );
        
        gsl_matrix_set(J, iter, 4, (1-difexp*(1-alpha+alpha/Z/Gamma*sinzg/powg) )/yi );
        gsl_matrix_set(J, iter, 5, 1/yi);

		//Punishment terms, to make constrains in parameter space. 
		if (alpha < 0)
		{
			gsl_matrix_set(J,iter,0, gsl_matrix_get(J, iter, 0)+ 2e5 *alpha);
		}

		if (alpha>1)
		{
			gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
		}

		if (D<0)
		{
			gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *D);
		}

		if (v<0)
		{
			gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *v);
		}

        if (Z<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *Z);
        }

        if (A<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif

#ifdef ISFRUNANDTUMBLE
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
	double* dataAry=((dataStruct *)sdata)->data;
	double* tau=((dataStruct *)sdata)->tau;
	double q=((dataStruct *)sdata)->q;

	double lambda=gsl_vector_get(para, 0);
	double v = gsl_vector_get(para, 1);
	double A = gsl_vector_get(para, 2);
	double B = gsl_vector_get(para, 3);

	double kv0=q*v;
	for (int iter = 0; iter<num_fit; ++iter)
	{
		int n=2;
		double incrementeven = 1;
		double result = incrementeven*gsl_sf_bessel_J0(kv0*tau[iter]);
		double bracket = lambda*lambda*tau[iter] / kv0;
		double incrementodd = sqrt(bracket);
		result += incrementodd*gsl_sf_bessel_j0(kv0*tau[iter])*sqrt(kv0*tau[iter]);
		double err=0;
        //Calculation of ISF by using the infinite series expansion, small trick is used for numerical stability.
		do
		{
			int n2 = n / 2;
			incrementeven *= bracket / (n - 1);
			result += incrementeven*gsl_sf_bessel_Jn(n2, kv0*tau[iter]);

			incrementodd *= bracket / n;
			err = incrementodd*gsl_sf_bessel_jl(n2, kv0*tau[iter])*sqrt(kv0*tau[iter]);
			result += err;
			n+=2;
		} while (abs(err/result)>precision);
		result *= exp(-lambda*tau[iter]);
		result = A*(1 - result) + B;

		gsl_vector_set(y, iter, result - dataAry[iter]);
    }
    return GSL_SUCCESS;
}

int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
	double* tau = ((dataStruct *)sdata)->tau;
	double q = ((dataStruct *)sdata)->q;

	double lambda = gsl_vector_get(para, 0);
	double v = gsl_vector_get(para, 1);
	double A = gsl_vector_get(para, 2);

	double kv0 = q*v;
	for (int iter = 0; iter<num_fit; ++iter)
	{
		double explt = exp(-lambda*tau[iter]);
		double bracket = lambda*lambda*tau[iter] / kv0;
		double relerr = 0;

		double incrementeven = 1;
		double j = gsl_sf_bessel_J0(kv0*tau[iter]);
		double f = j;
		double dfdl = tau[iter] * j;
		double dfdv = bracket*q*tau[iter]*j / 2.0; //j0 term

        //Calculation of Jacobian of ISF by using the infinite series expansion, small trick is used for numerical stability.

		double incrementodd = sqrt(bracket);
		j = incrementodd*gsl_sf_bessel_j0(kv0*tau[iter])*sqrt(kv0*tau[iter]);
		f += j;
		dfdl += (tau[iter] - 1.0 / lambda)*j;
		dfdv += j*q*tau[iter] / 4 * bracket;  //j1/2 term

		incrementeven *= bracket;
		j = gsl_sf_bessel_J1(kv0*tau[iter]);
		dfdv += q*tau[iter] * (bracket*bracket / 6.0 - 1)*j; //j1 term
		j *= incrementeven;
		f += j;
		dfdl += (tau[iter] - 2.0 / lambda)*j;

		incrementodd *= bracket / 2;
		j = incrementodd*gsl_sf_bessel_j1(kv0*tau[iter])*sqrt(kv0*tau[iter]);
		f += j;
		dfdl += j* (tau[iter] - 3.0 / lambda);
		dfdv += j*(lambda*lambda*tau[iter] * tau[iter] / 8.0 / v - q*q*v / lambda / lambda - 3.0 / 2.0 / v);

		dfdv -= lambda*q*tau[iter] * tau[iter] * gsl_sf_bessel_y0(kv0*tau[iter]); //y0 term
		int n = 4;
		do
		{
			int n2 = n / 2;
			incrementeven *= bracket / (n - 1);
			j = incrementeven*gsl_sf_bessel_Jn(n2, kv0*tau[iter]);
			f += j;
			dfdl += (tau[iter] - n / lambda)*j;
			dfdv += (lambda*lambda*tau[iter] * tau[iter] / 2.0/(n+1) / v - q*q*v*(n-1)/2.0 / lambda / lambda - n / 2.0 / v) *j;

			incrementodd *= bracket / n;
			j = incrementodd*gsl_sf_bessel_jl(n2, kv0*tau[iter])*sqrt(kv0*tau[iter]);
			f += j;
			relerr = abs(j / f);
			double ddfdl=j * (tau[iter] - (n + 1) / lambda);
			dfdl += ddfdl;
			if (ddfdl/dfdl>relerr)
			{
				relerr = ddfdl / dfdl;
			}
			double ddfdv = j* (lambda*lambda*tau[iter] * tau[iter] / 2.0 / (n + 2) / v - q*q*v*n / 2.0 / lambda / lambda - (n + 1) / 2.0 / v);
			dfdv += ddfdv;
			if (ddfdv / dfdv>relerr)
			{
				relerr = ddfdv / dfdv;
			}
			n += 2;
		} while (relerr>precision);
		f *= explt;
		dfdl *= explt;
		dfdv *= explt;

		gsl_matrix_set(J, iter, 0, A*dfdl);
		gsl_matrix_set(J, iter, 1, -A*dfdv);

		gsl_matrix_set(J, iter, 2, 1-f);
		gsl_matrix_set(J, iter, 3, 1);
    }
    return GSL_SUCCESS;
}

#endif

#ifdef ISFRUNANDTUMBLE_3D
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    double v0=gsl_vector_get(para, 0);
    double lambda=gsl_vector_get(para, 1);
    double A=gsl_vector_get(para, 2);
    double B=gsl_vector_get(para, 3);
    
    for (int iter = 0; iter<num_fit; ++iter)
    {
        double kv0=q*v0;
        double atanterm=atan(kv0/(tau[iter]+lambda));
        double yi=log(A*(1.0/tau[iter]-atanterm/(kv0-lambda*atanterm))+B/tau[iter]);

        double result = yi - dataAry[iter];
        //Punishment terms
        if (v0<0)
        {
            result += 1e5*v0*v0;
        }
        if (lambda<0)
        {
            result += 1e5*lambda*lambda;
        }
        if (A<0)
        {
            result += 1e5*A*A;
        }

        gsl_vector_set(y, iter, result);
    }

    return GSL_SUCCESS;
}

int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    double v0=gsl_vector_get(para, 0);
    double lambda=gsl_vector_get(para, 1);
    double A=gsl_vector_get(para, 2);
    double B=gsl_vector_get(para, 3);
    
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double kv0=q*v0;
        double lt=lambda+tau[iter];
        double atanterm=atan(kv0/lt);
        double lt2k2v02=lt*lt+kv0*kv0;
        double denominator2=kv0-lambda*atanterm;
        double denominator=lt2k2v02*denominator2*denominator2;
        double yi=A*(1.0/tau[iter]-atanterm/denominator2)+B/tau[iter];
        
        gsl_matrix_set(J, iter, 0, A*q*(lt2k2v02*atanterm-lt*kv0)/denominator/yi );
        
        gsl_matrix_set(J, iter, 1, A*(kv0*kv0-lt2k2v02*atanterm*atanterm)/denominator/yi );
        
        gsl_matrix_set(J, iter, 2, (1.0/tau[iter]-atanterm/denominator2)/yi );
        gsl_matrix_set(J, iter, 3, 1.0/tau[iter]/yi);

        //Punishment terms
        if (v0<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *v0);
        }
        if (lambda<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *lambda);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif

#ifdef ISFRunAndTumbleAndDiffusion
//The ISF is written to meet the API of GSL f function. sdata is the pointer to data structure defined by GSL. y is the return of the function.
int ISFfun(const gsl_vector* para, void* sdata, gsl_vector* y)
{
    double* dataAry=((dataStruct *)sdata)->data;
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    //Get the parameters.
    double alpha=gsl_vector_get(para, 0);
    double v0=gsl_vector_get(para, 1);
    double lambda=gsl_vector_get(para, 2);
    double D=gsl_vector_get(para, 3);
    double A=gsl_vector_get(para, 4);
    double B=gsl_vector_get(para, 5);
    
    double kv0=q*v0;
    double dq2=D*q*q;
    for (int iter = 0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double atanterm=atan(kv0/(tau[iter]+lambda+dq2));
        double yi=log(A*(1.0/tau[iter]-(1.0-alpha)/(dq2+tau[iter])-alpha*atanterm/(kv0-lambda*atanterm))+B/tau[iter]);

        double result = yi - dataAry[iter];
		//Punishment terms, to make constrains in parameter space. 
        if (v0<0)
        {
            result += 1e5*v0*v0;
        }
        if (lambda<0)
        {
            result += 1e5*lambda*lambda;
        }
        if (D<0)
        {
            result += 1e5*D*D;
        }
        if (alpha<0)
        {
            result += 1e5*alpha*alpha;
        }
        if (alpha>1)
        {
            result += 1e5*(alpha-1)*(alpha-1);
        }
        if (A<0)
        {
            result += 1e5*A*A;
        }

        gsl_vector_set(y, iter, result);
    }

    return GSL_SUCCESS;
}

//The function is written to meet the API of GSL df function. sdata is the pointer to data structure defined by GSL. J is the Jacobian, which is the return of the function.
int dISFfun(const gsl_vector* para, void* sdata, gsl_matrix* J)
{
    double* tau=((dataStruct *)sdata)->tau;
    double q=((dataStruct *)sdata)->q;
    
    double alpha=gsl_vector_get(para, 0);
    double v0=gsl_vector_get(para, 1);
    double lambda=gsl_vector_get(para, 2);
    double D=gsl_vector_get(para, 3);
    double A=gsl_vector_get(para, 4);
    double B=gsl_vector_get(para, 5);
    
    double kv0=q*v0;
    double dq2=D*q*q;
    for (int iter=0; iter<num_fit; ++iter)
    {
        //Temperary variables used for acceleration.
        double tdq2=tau[iter]+dq2;
        double lt=lambda+tdq2;
        double atanterm=atan(kv0/lt);
        double lt2k2v02=lt*lt+kv0*kv0;
        double denominator2=kv0-lambda*atanterm;
        double denominator=lt2k2v02*denominator2*denominator2;
        double yi=A*(1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*atanterm/denominator2)+B/tau[iter];
        
        gsl_matrix_set(J, iter, 0, A*(1.0/tdq2-atanterm/denominator2)/yi );
        gsl_matrix_set(J, iter, 1, A*q*alpha*(lt2k2v02*atanterm-lt*kv0)/denominator/yi );
        
        gsl_matrix_set(J, iter, 2, A*alpha*(kv0*kv0-lt2k2v02*atanterm*atanterm)/denominator/yi );

        gsl_matrix_set(J, iter, 3, A*q*q*((1.0-alpha)/tdq2/tdq2+alpha*kv0*kv0/denominator)/yi );

        gsl_matrix_set(J, iter, 4, (1.0/tau[iter]-(1.0-alpha)/tdq2-alpha*atanterm/denominator2)/yi );
        gsl_matrix_set(J, iter, 5, 1.0/tau[iter]/yi);

        //Punishment terms
        if (alpha<0)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *alpha);
        }
        if (alpha>1)
        {
            gsl_matrix_set(J, iter, 0, gsl_matrix_get(J, iter, 0) + 2e5 *(alpha-1));
        }
        if (v0<0)
        {
            gsl_matrix_set(J, iter, 1, gsl_matrix_get(J, iter, 1) + 2e5 *v0);
        }
        if (lambda<0)
        {
            gsl_matrix_set(J, iter, 2, gsl_matrix_get(J, iter, 2) + 2e5 *lambda);
        }
        if (D<0)
        {
            gsl_matrix_set(J, iter, 3, gsl_matrix_get(J, iter, 3) + 2e5 *D);
        }
        if (A<0)
        {
            gsl_matrix_set(J, iter, 4, gsl_matrix_get(J, iter, 4) + 2e5 *A);
        }
    }
    
    return GSL_SUCCESS;
}
#endif

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