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
    //ios::sync_with_stdio(false);
    stringstream arg;		//To read the argv
    arg << argv[1];
    
    ddm ddmExp;
    //    arg << "recover";
    
    if (arg.str() == "recover")		//recover from pre-calculated g(q,t) data
    {
        cout << "Loading datag.txt..." << '\n';
        ddmExp.recover();   //Loading g(q,t), q, and t data from files
    }
    else
    {
        string filePrefix = arg.str();
        //string filePrefix="cl-m1-ddm-01_X";
        
        cout << "Loading and calculating FFT..." << '\n';
        ddmExp.readAndFFT(filePrefix);
        
        cout << "Calculating average of square module for different tau... 0% finished." << '\n';
        ddmExp.averSqrModTau();
        ddmExp.cleanSeqk();     //Free memory

        cout << "Averaging on directions of q..." << '\n';
        ddmExp.aveQBilinear();
        ddmExp.cleankDiff();     //Free memory
       
        cout << "Printing unfitted data..." << '\n';
        ddmExp.printG();
    }
    
    cout << "Fitting..." << '\n';
    ddmExp.fitting();
    ddmExp.printFit();

    return 0;
}
