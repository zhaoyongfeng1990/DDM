//
//  print.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <fstream>
#include <iostream>

void ddm::printG()
{
    ofstream qfile("q.txt");
    for (int iter = 0; iter < qsize; ++iter)
    {
        qfile << qabs[iter] << endl;
    }
    qfile.close();
    ofstream datagfile("datag.txt");
    for (int iterq = 0; iterq < qsize; ++iterq)
    {
        for (int itertau = 0; itertau < num_fit; ++itertau)
        {
            datagfile << gsl_matrix_get(datag, iterq, itertau) << " ";
        }
        datagfile << endl;
    }
    datagfile.close();
}

void ddm::printFit()
{
    ofstream fitparafile("fitparafile.txt");
    ofstream fiterrfile("fiterrfile.txt");
    ofstream statusfile("status.txt");
    ofstream qstatusfile("statusq.txt");
    
    int cqsize=qsize-qIncreList[num_qCurve-1];
    int cnumOfPara=numOfPara+2*num_qCurve;
    
    for (int iterq=0; iterq<cqsize; ++iterq)
    {
        for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)
        {
            fitparafile << gsl_matrix_get(fittedPara, iterq, iterpara) << " ";
            fiterrfile << gsl_matrix_get(fitErr, iterq, iterpara) << " ";
        }
        fitparafile << endl;
        fiterrfile << endl;
        statusfile << iterq << ": q=" << qabs[iterq] << ", "<< gsl_strerror(status[iterq]) << endl;
        qstatusfile << status[iterq] << endl;
    }
    fitparafile.close();
    fiterrfile.close();
    statusfile.close();
    qstatusfile.close();
}