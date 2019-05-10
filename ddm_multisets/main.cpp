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
    
	ddm_Multisets ddmExp;
    //    arg << "recover";
    
    cout << "Loading datag.txt..." << '\n';
    ddmExp.recover();   //Loading g(q,t), q, and t data from files
    cout << "Fitting..." << '\n';
    ddmExp.fitting();
    ddmExp.printFit();

    return 0;
}
