//
//  loadFiles.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <sstream>
#include <fstream>
#include <omp.h>

void readTiff(const string tifName, gsl_matrix* data, const int dimx, const int dimy);
//Reading tiff images.
void readSim(const string fileName, gsl_matrix* data, const int dimx, const int dimy);
//Read simulation data from ASCII file.

//The data structures for reading tiff
typedef struct tagTIFFILEHEAD{
    unsigned short biOrder;
    unsigned short version;
    uint32_t offsetIFD;
} TIFFILEHEAD;

typedef struct tagTIFFDE{
    unsigned short tag;
    unsigned short type;
    uint32_t length;
    uint32_t valueOffset;
} TIFFDE;

typedef struct tagTIFFIFD{
    unsigned short DEC;
    TIFFDE* DE;
    uint32_t offset;
    tagTIFFIFD* nextIFD;
} TIFFIFD;

void ddm::readAndFFT(const string filePrefix)
{
    int cdimy=dimy;
    int cdimx=dimx;
    int cdimkx=dimkx;
#pragma omp parallel
    {
        gsl_matrix* fftMatrix = gsl_matrix_alloc(cdimy, cdimx);
        gsl_matrix_complex* resultMatrix = gsl_matrix_complex_alloc(cdimy, cdimkx);
        //dimk=dim/2+1, FFTW only return around half points for r2c FFT to save memory usage.
        
        fftw_plan fft2plan;		//FFT need to make a plan before reading data.
#pragma omp critical (make_plan)
        fft2plan = fftw_plan_dft_r2c_2d(cdimy, cdimx, fftMatrix->data, (fftw_complex *)resultMatrix->data, FFTW_MEASURE);
        
#pragma omp for
        for (int iter = 0; iter < numOfSeq; ++iter)
        {
            imageSeqk[iter] = gsl_matrix_complex_calloc(cdimy, cdimkx);
            if (filePrefix=="simulation")
            {
                stringstream fileName;
                fileName << iter << ".txt";
                readSim(fileName.str(), fftMatrix, cdimx, cdimy);
            }
            else
            {
                stringstream fileName;
                fileName << filePrefix << iter+1 << ".tif";		//Making path name
                readTiff(fileName.str(), fftMatrix, cdimx, cdimy);
            }
            fftw_execute(fft2plan);
            
            gsl_matrix_complex_memcpy(imageSeqk[iter], resultMatrix);
        }
        gsl_matrix_free(fftMatrix);
        gsl_matrix_complex_free(resultMatrix);
        fftw_destroy_plan(fft2plan);
    }
}

void ddm::recover()
{
    qabs.reserve((dimkx - 1)*(dimky - 1)/2);	//the number of effective q will not exceed (dimk - 1)*(dimk - 1)/2.
    ifstream qfile("q.txt");
    while (!qfile.eof())
    {
        double tempq;
        qfile >> tempq;
        qabs.push_back(tempq);
    }
    
    qabs.pop_back();
    qsize=(int)qabs.size();
    qfile.close();
    
    tau.reserve(numOfDiff);
    ifstream taufile("tau.txt");
    while (!taufile.eof())
    {
        double temptau;
        taufile >> temptau;
        tau.push_back(temptau);
    }
    tau.pop_back();
    num_fit=(int)tau.size();
    taufile.close();
    
    datag = gsl_matrix_alloc(qsize, num_fit);		//Allocate memory for g(q,t) matrix.
    ifstream datagfile("datag.txt");
    for (int iterq = 0; iterq < qsize; ++iterq)		//Reading g(q,t) matrix.
    {
        for (int itertau = 0; itertau < num_fit; ++itertau)
        {
            double tempdata;
            datagfile >> tempdata;
            gsl_matrix_set(datag, iterq, itertau, tempdata);
        }
    }
    datagfile.close();
}

//Reading tif files
void readTiff(const string tifName, gsl_matrix* data, const int dimx, const int dimy)
{
    int width=dimx;	//Width and height of the image equals to dim (e.g. =512 pixel)
    int height=dimy;
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
}

void readSim(const string fileName, gsl_matrix* data, const int dimx, const int dimy)	//Read simulation data from ASCII file.
{
    ifstream infile(fileName.c_str());
    for (int iterx=0; iterx<dimy; ++iterx)
    {
        for (int itery=0; itery<dimx; ++itery)
        {
            double temp;
            infile >> temp;
            gsl_matrix_set(data, iterx, itery, temp);
        }
    }
}