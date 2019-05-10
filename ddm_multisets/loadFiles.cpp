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
#include <iostream>
#include <omp.h>

void ddm_Multisets::recover()
{
    //Loading table of q
    ifstream qfile("q_low.txt");
    while (!qfile.eof())
    {
        double tempq;
        qfile >> tempq;
        qabs_low.push_back(tempq);
    }
    //The last line must be '\n', drop it. (So there must be a '\n' at the end of the file! )
    qabs_low.pop_back();
    qsize_low=(int)qabs_low.size();
    qfile.close();

	qfile.open("q_high.txt");
	while (!qfile.eof())
	{
		double tempq;
		qfile >> tempq;
		qabs_high.push_back(tempq);
	}
	//The last line must be '\n', drop it. (So there must be a '\n' at the end of the file! )
	qabs_high.pop_back();
	qsize_high = (int)qabs_high.size();
	qfile.close();

    //Loading table of tau
    ifstream taufile("tau_low.txt");
    while (!taufile.eof())
    {
        double temptau;
        taufile >> temptau;
        tau_low.push_back(temptau);
    }
    //The last line must be '\n', drop it. (So there must be a '\n' at the end of the file! )
    tau_low.pop_back();
    num_fit_low=(int)tau_low.size();
    taufile.close();

	taufile.open("tau_high.txt");
	while (!taufile.eof())
	{
		double temptau;
		taufile >> temptau;
		tau_high.push_back(temptau);
	}
	//The last line must be '\n', drop it. (So there must be a '\n' at the end of the file! )
	tau_high.pop_back();
	num_fit_high = (int)tau_high.size();
	taufile.close();
    
    datag_low = gsl_matrix_alloc(qsize_low, num_fit_low);		//Allocate memory for g(q,t) matrix.
	datag_high = gsl_matrix_alloc(qsize_high, num_fit_high);
	ifstream datagfile("datag_low.txt");
    for (int iterq = 0; iterq < qsize_low; ++iterq)		//Reading g(q,t) matrix.
    {
        for (int itertau = 0; itertau < num_fit_low; ++itertau)
        {
            double tempdata;
            datagfile >> tempdata;
            gsl_matrix_set(datag_low, iterq, itertau, tempdata);
        }
    }
    datagfile.close();

	datagfile.open("datag_high.txt");
	for (int iterq = 0; iterq < qsize_high; ++iterq)		//Reading g(q,t) matrix.
	{
		for (int itertau = 0; itertau < num_fit_high; ++itertau)
		{
			double tempdata;
			datagfile >> tempdata;
			gsl_matrix_set(datag_high, iterq, itertau, tempdata);
		}
	}
	datagfile.close();
}
