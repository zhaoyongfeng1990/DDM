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

#define MIN(a,b) a>b?b:a
//Print unfitted data

//Print fitting result
void ddm_Multisets::printFit()
{
    ofstream fitparafile("fitparafile.txt");
    ofstream fiterrfile("fiterrfile.txt");
    ofstream statusfile("qcurves.txt");
    ofstream qstatusfile("statusq.txt");

	const int cqsize = MIN(qsize_low - qIncreList_low[num_qCurve_low - 1], qsize_high - qIncreList_high[num_qCurve_high - 1]);
    int cnumOfPara=numOfPara+2*(num_qCurve_low+num_qCurve_high);
    
    for (int iterq=0; iterq<cqsize; ++iterq)
    {
        for (int iterpara=0; iterpara<cnumOfPara; ++iterpara)
        {
            fitparafile << gsl_matrix_get(fittedPara, iterq, iterpara) << " ";
            fiterrfile << gsl_matrix_get(fitErr, iterq, iterpara) << " ";
        }
        fitparafile << '\n';
        fiterrfile << '\n';
        qstatusfile << status[iterq] << '\n';
    }
	for (int iterq = 0; iterq < cqsize; ++iterq)
	{
		for (int iterqc = 0; iterqc < num_qCurve_low; ++iterqc)
		{
			statusfile << qabs_low[iterq + qIncreList_low[iterqc]] << " ";
		}
		for (int iterqc = 0; iterqc < num_qCurve_high; ++iterqc)
		{
			statusfile << qabs_high[iterq + qIncreList_high[iterqc]] << " ";
		}
		statusfile << "\n";
	}
    fitparafile.close();
    fiterrfile.close();
    statusfile.close();
    qstatusfile.close();
}

//Print arbitrary gsl_matrix m for debug
void ddm_Multisets::printdebugM(gsl_matrix* m, const string filename)
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