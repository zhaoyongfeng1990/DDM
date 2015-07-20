//
//  aveQ.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <iostream>

//int find(const vector<int>& vec, const int value);
////Find the position of a particular member in a vector.
//int quickFind(const vector<int>& vec, const int value);
////Find the position of a particular member in a sorted vector with increasing order, using dichotomy.

//Average all directions of q by bilinear interpolation. This is done by integrating the bilinear interpolation function in every lattice, and then dividing by the length of arc.
void ddm::aveQBilinear()
{
    double cqmin=qmin;
    double cqmax=qmax;
    double cqstep=qstep;
    double cdqx=dqx;
    double cdqy=dqy;
    int cnum_fit=num_fit;
    int cdimy=dimy;
    int cdimkx=dimkx;
    int cdimky=dimky;
    double cqsize=ceil((cqmax-cqmin)/cqstep); //qsize is the number of different q value samples.
    qabs.resize(cqsize);     //qabs is the absolute value of q
    for (int iter=0; iter<cqsize; ++iter)
    {
        qabs[iter]=iter*cqstep+cqmin;    //For interpolation.
    }
    datag = gsl_matrix_calloc(cqsize, cnum_fit);
    gsl_matrix* count = gsl_matrix_alloc(cqsize, cnum_fit);
    gsl_matrix_set_zero(datag);
    gsl_matrix_set_zero(count);     //Number of elements
    
#pragma omp parallel for
    for (int itertau = 0; itertau < cnum_fit; ++itertau)
    {
        for (int iterrow = 1; iterrow < cdimky-1; ++iterrow)
        {
            int kx1 = iterrow;
            int kx2 = iterrow+1;
            int refkx1=(cdimy-kx1)%cdimy;
            int refkx2=(cdimy-kx2)%cdimy;
            for (int itercol = 1; itercol < cdimkx-1; ++itercol)
            {
                int ky1 = itercol;
                int ky2 = itercol+1;
                
                double dist1=sqrt(kx1*kx1*cdqy*cdqy+ky1*ky1*cdqx*cdqx)/cqstep;
                double dist2=sqrt(kx2*kx2*cdqy*cdqy+ky2*ky2*cdqx*cdqx)/cqstep;
                
                int maxqidx=ceil(dist2-cqmin/cqstep);
                int minqidx=ceil(dist1-cqmin/cqstep);
                
                double dist3=sqrt(kx2*kx2*cdqy*cdqy+ky1*ky1*cdqx*cdqx)/cqstep;
                double dist4=sqrt(kx1*kx1*cdqy*cdqy+ky2*ky2*cdqx*cdqx)/cqstep;
                
                for (int iterq=minqidx; iterq<maxqidx; ++iterq)
                {
                    double px[2];
                    double py[2];
                    
                    double cq=iterq+cqmin/cqstep;
                    
                    //cout << cq << '\n';
                    
                    if(dist3>=cq)
                    {
                        py[0]=ky1*cdqx/cqstep;
                        px[0]=sqrt(cq*cq-py[0]*py[0]);
                    }
                    else
                    {
                        px[0]=kx2*cdqy/cqstep;
                        py[0]=sqrt(cq*cq-px[0]*px[0]);
                    }
                    
                    if (dist4>=cq)
                    {
                        px[1]=kx1*cdqy/cqstep;
                        py[1]=sqrt(cq*cq-px[1]*px[1]);
                    }
                    else
                    {
                        py[1]=ky2*cdqx/cqstep;
                        px[1]=sqrt(cq*cq-py[1]*py[1]);
                    }
                    
                    double dist=(px[0]-px[1])*(px[0]-px[1])+(py[0]-py[1])*(py[0]-py[1]);
                    double cosdt=1-dist/2/cq/cq;
                    double dtheta=acos(cosdt);
                    
                    double cost1=px[0]/cq;
                    double sint1=py[0]/cq;
                    double cost2=px[1]/cq;
                    double sint2=py[1]/cq;
                    double dcost=cost2-cost1;
                    double dsint=sint2-sint1;
                    double dcos2t=2*(cost2*cost2-cost1*cost1);
                    
                    double u11=gsl_matrix_get(imagekDiff[itertau], kx1, ky1);
                    double u12=gsl_matrix_get(imagekDiff[itertau], kx1, ky2);
                    double u21=gsl_matrix_get(imagekDiff[itertau], kx2, ky1);
                    double u22=gsl_matrix_get(imagekDiff[itertau], kx2, ky2);
                    
                    cq*=cqstep;
                    double arc=gsl_matrix_get(datag, iterq, itertau);
                    arc+=u11*dtheta-(cq*dcost+ky1*cdqx*dtheta)*(u12-u11)/cdqx+(cq*dsint-kx1*cdqy*dtheta)*(u21-u11)/cdqy+(kx1*cdqy*cq*dcost+kx1*ky1*cdqx*cdqy*dtheta-cq*cq*dcos2t/4-ky1*cdqx*cq*dsint)/cdqx/cdqy*(u22-u21-u12+u11);
                    
                    //cout << arc << '\n';
                    
                    u11=gsl_matrix_get(imagekDiff[itertau], refkx1, ky1);
                    u12=gsl_matrix_get(imagekDiff[itertau], refkx1, ky2);
                    u21=gsl_matrix_get(imagekDiff[itertau], refkx2, ky1);
                    u22=gsl_matrix_get(imagekDiff[itertau], refkx2, ky2);
                    
                    arc+=u11*dtheta-(cq*dcost+ky1*cdqx*dtheta)*(u12-u11)/cdqx+(cq*dsint-kx1*cdqy*dtheta)*(u21-u11)/cdqy+(kx1*cdqy*cq*dcost+kx1*ky1*cdqx*cdqy*dtheta-cq*cq*dcos2t/4-ky1*cdqx*cq*dsint)/cdqx/cdqy*(u22-u21-u12+u11);
                    
                    gsl_matrix_set(datag, iterq, itertau, arc);
                    double angle=gsl_matrix_get(count, iterq, itertau);
                    angle+=dtheta*2;
                    gsl_matrix_set(count, iterq, itertau, angle);
                }
            }
        }
    }
    gsl_matrix_div_elements(datag, count);		//Average
    qsize=cqsize;
}

void ddm::aveQBicubic()
{
    qsize=ceil(qmax/qstep); //qsize is the number of different q value samples.
    qabs.resize(qsize);     //qabs is the absolute value of q
    for (int iter=0; iter<qsize; ++iter)
    {
        qabs[iter]=iter*qstep;    //For interpolation.
    }
    datag = gsl_matrix_calloc(qsize, numOfDiff);
    gsl_matrix* count = gsl_matrix_alloc(qsize, numOfDiff);
    gsl_matrix_set_zero(datag);
    gsl_matrix_set_zero(count);     //Number of elements
    
#pragma omp parallel for
    for (int itertau = 0; itertau < numOfDiff; ++itertau)
    {
        for (int iterrow = 0; iterrow < dimky-1; ++iterrow)
        {
            int kx1 = iterrow;
            int kx2 = iterrow+1;
            int refkx1=(dimy-kx1)%dimy;
            int refkx2=(dimy-kx2)%dimy;
            for (int itercol = 0; itercol < dimkx-1; ++itercol)
            {
                int ky1 = itercol;
                int ky2 = itercol+1;
                
                double dist1=sqrt(kx1*kx1*dqy*dqy+ky1*ky1*dqx*dqx)/qstep;
                double dist2=sqrt(kx2*kx2*dqy*dqy+ky2*ky2*dqx*dqx)/qstep;
                
                int maxqidx=ceil(dist2);
                int minqidx=ceil(dist1);
                
                double dist3=sqrt(kx2*kx2*dqy*dqy+ky1*ky1*dqx*dqx)/qstep;
                double dist4=sqrt(kx1*kx1*dqy*dqy+ky2*ky2*dqx*dqx)/qstep;
                
                for (int iterq=minqidx; iterq<maxqidx; ++iterq)
                {
                    double px[2];
                    double py[2];
                    
                    if(dist3>=iterq)
                    {
                        py[0]=ky1*dqx/qstep;
                        px[0]=sqrt(iterq*iterq-py[0]*py[0]);
                    }
                    else
                    {
                        px[0]=kx2*dqy/qstep;
                        py[0]=sqrt(iterq*iterq-px[0]*px[0]);
                    }
                    
                    if (dist4>=iterq)
                    {
                        px[1]=kx1*dqy/qstep;
                        py[1]=sqrt(iterq*iterq-px[1]*px[1]);
                    }
                    else
                    {
                        py[1]=ky2*dqx/qstep;
                        px[1]=sqrt(iterq*iterq-py[1]*py[1]);
                    }
                    
                    double dist=(px[0]-px[1])*(px[0]-px[1])+(py[0]-py[1])*(py[0]-py[1]);
                    double cosdt=1-dist/2/iterq/iterq;
                    double dtheta=acos(cosdt);
                    
                    double cost1=px[0]/iterq;
                    double sint1=py[0]/iterq;
                    double cost2=px[1]/iterq;
                    double sint2=py[1]/iterq;
                    double dcost=cost2-cost1;
                    double dsint=sint2-sint1;
                    double dcos2t=2*(cost2*cost2-cost1*cost1);
                    
                    double u11=gsl_matrix_get(imagekDiff[itertau], kx1, ky1);
                    double u12=gsl_matrix_get(imagekDiff[itertau], kx1, ky2);
                    double u21=gsl_matrix_get(imagekDiff[itertau], kx2, ky1);
                    double u22=gsl_matrix_get(imagekDiff[itertau], kx2, ky2);
                    
                    double arc=gsl_matrix_get(datag, iterq, itertau);
                    arc+=u11*dtheta-(qabs[iterq]*dcost+ky1*dqx*dtheta)*(u12-u11)/dqx+(qabs[iterq]*dsint-kx1*dqy*dtheta)*(u21-u11)/dqy+(kx1*dqy*qabs[iterq]*dcost+kx1*ky1*dqx*dqy*dtheta-qabs[iterq]*qabs[iterq]*dcos2t/4-ky1*dqx*qabs[iterq]*dsint)/dqx/dqy*(u22-u21-u12+u11);
                    
                    u11=gsl_matrix_get(imagekDiff[itertau], refkx1, ky1);
                    u12=gsl_matrix_get(imagekDiff[itertau], refkx1, ky2);
                    u21=gsl_matrix_get(imagekDiff[itertau], refkx2, ky1);
                    u22=gsl_matrix_get(imagekDiff[itertau], refkx2, ky2);
                    
                    arc+=u11*dtheta-(qabs[iterq]*dcost+ky1*dqx*dtheta)*(u12-u11)/dqx+(qabs[iterq]*dsint-kx1*dqy*dtheta)*(u21-u11)/dqy+(kx1*dqy*qabs[iterq]*dcost+kx1*ky1*dqx*dqy*dtheta-qabs[iterq]*qabs[iterq]*dcos2t/4-ky1*dqx*qabs[iterq]*dsint)/dqx/dqy*(u22-u21-u12+u11);
                    
                    gsl_matrix_set(datag, iterq, itertau, arc);
                    double angle=gsl_matrix_get(count, iterq, itertau);
                    angle+=dtheta*2;
                    gsl_matrix_set(count, iterq, itertau, angle);
                }
            }
        }
        gsl_matrix_set(datag, 0, itertau, gsl_matrix_get(imagekDiff[itertau], 0, 0));
        gsl_matrix_set(count, 0, itertau, 1);
    }
    gsl_matrix_div_elements(datag, count);		//Average
}
