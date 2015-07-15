//
//  aveQ.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"

//int find(const vector<int>& vec, const int value);
////Find the position of a particular member in a vector.
//int quickFind(const vector<int>& vec, const int value);
////Find the position of a particular member in a sorted vector with increasing order, using dichotomy.

//Average all directions of q by bilinear interpolation. This is done by integrating the bilinear interpolation function in every lattice, and then dividing by the length of arc.
void ddm::aveQBilinear()
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
        for (int iterrow = 1; iterrow < dimky-1; ++iterrow)
        {
            int kx1 = iterrow;
            int kx2 = iterrow+1;
            int refkx1=(dimy-kx1)%dimy;
            int refkx2=(dimy-kx2)%dimy;
            for (int itercol = 1; itercol < dimkx-1; ++itercol)
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
                    double dt=acos(cosdt);
                    
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
                    arc+=u11*dt-(qabs[iterq]*dcost+ky1*dqx*dt)*(u12-u11)/dqx+(qabs[iterq]*dsint-kx1*dqy*dt)*(u21-u11)/dqy+(kx1*dqy*qabs[iterq]*dcost+kx1*ky1*dqx*dqy*dt-qabs[iterq]*qabs[iterq]*dcos2t/4-ky1*dqx*qabs[iterq]*dsint)/dqx/dqy*(u22-u21-u12+u11);
                    
                    u11=gsl_matrix_get(imagekDiff[itertau], refkx1, ky1);
                    u12=gsl_matrix_get(imagekDiff[itertau], refkx1, ky2);
                    u21=gsl_matrix_get(imagekDiff[itertau], refkx2, ky1);
                    u22=gsl_matrix_get(imagekDiff[itertau], refkx2, ky2);
                    
                    arc+=u11*dt-(qabs[iterq]*dcost+ky1*dqx*dt)*(u12-u11)/dqx+(qabs[iterq]*dsint-kx1*dqy*dt)*(u21-u11)/dqy+(kx1*dqy*qabs[iterq]*dcost+kx1*ky1*dqx*dqy*dt-qabs[iterq]*qabs[iterq]*dcos2t/4-ky1*dqx*qabs[iterq]*dsint)/dqx/dqy*(u22-u21-u12+u11);
                    
                    gsl_matrix_set(datag, iterq, itertau, arc);
                    double angle=gsl_matrix_get(count, iterq, itertau);
                    angle+=dt*2;
                    gsl_matrix_set(count, iterq, itertau, angle);
                }
            }
        }
        gsl_matrix_set(datag, 0, itertau, gsl_matrix_get(imagekDiff[itertau], 0, 0));
        gsl_matrix_set(count, 0, itertau, 1);
    }
    gsl_matrix_div_elements(datag, count);		//Average
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
                    double dt=acos(cosdt);
                    
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
                    arc+=u11*dt-(qabs[iterq]*dcost+ky1*dqx*dt)*(u12-u11)/dqx+(qabs[iterq]*dsint-kx1*dqy*dt)*(u21-u11)/dqy+(kx1*dqy*qabs[iterq]*dcost+kx1*ky1*dqx*dqy*dt-qabs[iterq]*qabs[iterq]*dcos2t/4-ky1*dqx*qabs[iterq]*dsint)/dqx/dqy*(u22-u21-u12+u11);
                    
                    u11=gsl_matrix_get(imagekDiff[itertau], refkx1, ky1);
                    u12=gsl_matrix_get(imagekDiff[itertau], refkx1, ky2);
                    u21=gsl_matrix_get(imagekDiff[itertau], refkx2, ky1);
                    u22=gsl_matrix_get(imagekDiff[itertau], refkx2, ky2);
                    
                    arc+=u11*dt-(qabs[iterq]*dcost+ky1*dqx*dt)*(u12-u11)/dqx+(qabs[iterq]*dsint-kx1*dqy*dt)*(u21-u11)/dqy+(kx1*dqy*qabs[iterq]*dcost+kx1*ky1*dqx*dqy*dt-qabs[iterq]*qabs[iterq]*dcos2t/4-ky1*dqx*qabs[iterq]*dsint)/dqx/dqy*(u22-u21-u12+u11);
                    
                    gsl_matrix_set(datag, iterq, itertau, arc);
                    double angle=gsl_matrix_get(count, iterq, itertau);
                    angle+=dt*2;
                    gsl_matrix_set(count, iterq, itertau, angle);
                }
            }
        }
        gsl_matrix_set(datag, 0, itertau, gsl_matrix_get(imagekDiff[itertau], 0, 0));
        gsl_matrix_set(count, 0, itertau, 1);
    }
    gsl_matrix_div_elements(datag, count);		//Average
}

void ddm::aveQSort()
{
    //    vector<int> q;
    //    q.reserve((dimk - 1)*(dimk - 1) / 2);
    //    q.push_back(0);
    //    for (int iteri = 1; iteri < dimk; ++iteri)
    //    {
    //        for (int iterj = 0; iterj <= iteri; ++iterj) {
    //            int index = iteri*iteri + iterj*iterj;
    //            if (find(q, index) == -1)
    //            {
    //                q.push_back(index);
    //            }
    //        }
    //    }
    //
    //    qsize = (int)q.size();
    //    //sorting q
    //    int stackIter = 0;
    //    int* stack = new int[qsize];
    //    stack[stackIter] = 0;
    //    ++stackIter;
    //    stack[stackIter] = qsize - 1;
    //
    //    while (stackIter > -1)
    //    {
    //        int right = stack[stackIter];
    //        --stackIter;
    //        int left = stack[stackIter];
    //        --stackIter;
    //        int empty = left - 1;
    //        double pivotvalue = q[right];
    //        for (int64_t iter = left; iter < right; ++iter)
    //        {
    //            if (q[iter] < pivotvalue)
    //            {
    //                ++empty;
    //                int temp = q[iter];
    //                q[iter] = q[empty];
    //                q[empty] = temp;
    //            }
    //        }
    //        int temp = q[right];
    //        q[right] = q[empty + 1];
    //        q[empty + 1] = temp;
    //        if (left < empty)
    //        {
    //            ++stackIter;
    //            stack[stackIter] = left;
    //            ++stackIter;
    //            stack[stackIter] = empty;
    //        }
    //        if (empty + 2 < right)
    //        {
    //            ++stackIter;
    //            stack[stackIter] = empty + 2;
    //            ++stackIter;
    //            stack[stackIter] = right;
    //        }
    //    }
    //    delete[] stack;
    //
    //    datag = gsl_matrix_calloc(qsize, numOfDiff);
    //    gsl_matrix* count = gsl_matrix_alloc(qsize, numOfDiff);
    //    gsl_matrix_set_zero(datag);
    //    gsl_matrix_set_zero(count);
    //#pragma omp parallel for
    //    for (int itertau = 0; itertau < numOfDiff; ++itertau)
    //    {
    //        for (int iterrow = 0; iterrow < dim; ++iterrow)
    //        {
    //            int kx;
    //            if (iterrow < dimk)
    //                kx = iterrow;
    //            else
    //                kx = dim - iterrow;
    //            for (int itercol = 0; itercol < dimk; ++itercol)
    //            {
    //                int ky = itercol;
    //                int factor = 2;
    //                if (ky % (dimk - 1) == 0)
    //                    factor = 1;
    //                int k2 = kx*kx + ky*ky;
    //                int index = quickFind(q, k2);
    //
    //                double temp = gsl_matrix_get(imagekDiff[itertau], iterrow, itercol)*factor;
    //                temp += gsl_matrix_get(datag, index, itertau);
    //                gsl_matrix_set(datag, index, itertau, temp);
    //                double tempindex = gsl_matrix_get(count, index, itertau) + factor;
    //                gsl_matrix_set(count, index, itertau, tempindex);
    //            }
    //        }
    //    }
    //#pragma omp parallel for
    //    for (int iter = 0; iter < numOfDiff; ++iter)
    //        gsl_matrix_free(imagekDiff[iter]);
    //    imagekDiff.clear();
}

void ddm::aveQConcentric()
{
//    qsize=ceil(qmax/qstep);
//    qabs.resize(qsize);     //qabs is the absolute value of q
//    for (int iter=0; iter<qsize; ++iter)
//    {
//        //qabs[iter]=(iter+0.5)*qstep;    //For averaging. The midpoint of the cirque
//        qabs[iter]=iter*qstep;    //For interpolation.
//    }
//    
//    datag = gsl_matrix_calloc(qsize, numOfDiff);
//    gsl_matrix* count = gsl_matrix_alloc(qsize, numOfDiff);
//    gsl_matrix_set_zero(datag);
//    gsl_matrix_set_zero(count);     //Number of elements

//#pragma omp parallel for
//    for (int itertau = 0; itertau < numOfDiff; ++itertau)
//    {
//        for (int iterrow = 0; iterrow < dim; ++iterrow)
//        {
//            int kx;
//            if (iterrow < dimk)		//kx=[0, 1, ..., dim/2, 1-dim/2, 2-dim/2, ..., -1]
//                kx = iterrow;
//            else
//                kx = iterrow - dim;
//            for (int itercol = 0; itercol < dimk; ++itercol)
//            {
//                int ky = itercol;
//                int factor = 2;     //If the data is not the q=0 or N-1, it appears twice in q-space for a_i=a_{dim-i}.
//                if (ky % (dimk - 1) == 0)
//                    factor = 1;
//                
//                double k2 = kx*kx + ky*ky;
//                int index = floor(sqrt(k2)*dq/qstep);
//                
//                //This is fast enough. No need to speed up.
//                double temp = gsl_matrix_get(imagekDiff[itertau], iterrow, itercol)*factor;
//                temp += gsl_matrix_get(datag, index, itertau);
//                gsl_matrix_set(datag, index, itertau, temp);		//datag[index, itertau]+=imagekDiff[itertau][tierrow, itercol]*factor
//                double tempcount = gsl_matrix_get(count, index, itertau) + factor;
//                gsl_matrix_set(count, index, itertau, tempcount);	//count[index, itertau]+=factor
//            }
//        }
//    }
//    gsl_matrix_div_elements(datag, count);		//Average
//#pragma omp parallel for
//    for (int iter = 0; iter < numOfDiff; ++iter)
//        gsl_matrix_free(imagekDiff[iter]);
//    imagekDiff.clear();
}

///////////////////////////
////Functions not in use.//
///////////////////////////
//
////Find the position of a particular member in a sorted vector with increasing order, using dichotomy.
//int quickFind(const vector<int>& vec, const int value)
//{
//    int left=0;
//    int right=(int)vec.size();
//    while(right-left>1)
//    {
//        int middle=(left+right)/2;
//        if (vec[middle]==value)
//        {
//            return middle;
//        }
//        else if (vec[middle]>value)
//        {
//            right=middle;
//        }
//        else
//        {
//            left=middle;
//        }
//    }
//    if (vec[left]==value)
//        return left;
//    else if (vec[right]==value)
//        return right;
//    else
//        return -1;
//}
//
////Find the position of a particular member in a vector.
//int find(const vector<int>& vec, const int value)
//{
//    for (int iter=0; iter<vec.size(); ++iter)
//    {
//        if (vec[iter]==value)
//            return iter;
//    }
//    return -1;
//}