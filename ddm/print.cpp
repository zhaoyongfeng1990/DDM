//
//  print.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 15/6/22.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <fstream>
//#include <iostream>

//Print unfitted data
void ddm::printG()
{
    ofstream qfile("q.txt");
    for (int iter = 0; iter < qsize; ++iter)
    {
        qfile << qabs[iter] << '\n';
    }
    qfile.close();
    ofstream datagfile("datag.txt");
    for (int iterq = 0; iterq < qsize; ++iterq)
    {
        for (int itertau = 0; itertau < num_fit; ++itertau)
        {
            datagfile << gsl_matrix_get(datag, iterq, itertau) << " ";
        }
        datagfile << '\n';
    }
    datagfile.close();
    ofstream taufile("tau.txt");
    for (int itertau = 0; itertau < num_fit; ++itertau)
    {
        taufile << tau[itertau] << '\n';
    }
}

//Print fitting result
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
        fitparafile << '\n';
        fiterrfile << '\n';
        statusfile << iterq << ": q=" << qabs[iterq] << ", "<< gsl_strerror(status[iterq]) << '\n';
        qstatusfile << status[iterq] << '\n';
    }
    fitparafile.close();
    fiterrfile.close();
    statusfile.close();
    qstatusfile.close();
}

//Print arbitrary gsl_matrix m for debug
void ddm::printdebugM(gsl_matrix* m, const string filename)
{
    ofstream outfile(filename);
    for (int iter=0; iter<m->size1; ++iter)
    {
        for (int iterc=0; iterc<m->size2; ++iterc)
        {
            outfile << gsl_matrix_get(m, iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}