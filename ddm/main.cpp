//
//  main.cpp
//  ddm
//
//  Created by Zhao Yongfeng on 14/11/3.
//  Copyright (c) 2014年 ZYF. All rights reserved.
//

#include "ddm.h"
#include <sstream>
#include <iostream>
#include <omp.h>

//The name prefix of the images should be given in bash as argv. (e.g.: IM_340_01_X)
//Or, give "simulation" can deal with the simulation data.
//Or, give "recover" can do the fit by reading datag.txt and q.txt directly, save doing time-consuming FFT and averaging.
int main(int argc, const char * argv[])
{
    omp_set_num_threads(OMP_NUM_THREADS);     //Number of threads
    stringstream arg;		//To read the argv
    arg << argv[1];
    
    ddm ddmExp;
    //    arg << "recover";
    
    if (arg.str() == "recover")		//recover from pre-calculated g(q,t) data
    {
        cout << "Loading datag.txt..." << endl;
        ddmExp.recover();
    }
    else
    {
        string filePrefix = arg.str();
        //string filePrefix="cl-m1-ddm-01_X";
        
        cout << "Loading and calculating FFT..." << endl;
        ddmExp.readAndFFT(filePrefix);
        
        cout << "Calculating average of square module for different tau... 0% finished." << endl;
        ddmExp.averSqrModTau();
        ddmExp.cleanSeqk();

        cout << "Averaging on directions of q..." << endl;
        ddmExp.aveQBilinear();
        ddmExp.cleankDiff();
       
        cout << "Printing unfitted data..." << endl;
        ddmExp.printG();
    }
    
#ifdef NeedLaplaceTrans
    cout << "Numerical Laplace transforming..." << endl;
    ddmExp.LaplaceTrans();
    ddmExp.printGs();
#endif
    
    cout << "Fitting..." << endl;
    ddmExp.fitting_estRange();
    //ddmExp.printFit();

    return 0;
}
