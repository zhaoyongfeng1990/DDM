#!/bin/bash

#  config.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

omp_num=8
nqcurve=1

qincre=

dimx=512
dimy=512
numOfSeq=4500
numOfDiff=4000
dx=6.5/4

qmin=0.01
qstep=0.01
dt=0.01
maxIter=200
alphaGuess=0.4
DGuess=0.2
vbarGuess=6
lambdaGuess=0.5
ZGuess=15
sigmaGuess=3

if [ $nqcurve == 1 ]
then
    paraList="$omp_num""\n""$nqcurve""\n""$dimx""\n""$dimy""\n""$numOfSeq""\n""$numOfDiff""\n""$dx""\n""$qmin""\n""$qstep""\n""$dt""\n""$maxIter""\n""$alphaGuess""\n""$DGuess""\n""$vbarGuess""\n""$lambdaGuess""\n""$ZGuess""\n""$sigmaGuess""\n"
else
    paraList="$omp_num""\n""$nqcurve""\n""$qincre""\n""$dimx""\n""$dimy""\n""$numOfSeq""\n""$numOfDiff""\n""$dx""\n""$qmin""\n""$qstep""\n""$dt""\n""$maxIter""\n""$alphaGuess""\n""$DGuess""\n""$vbarGuess""\n""$lambdaGuess""\n""$ZGuess""\n""$sigmaGuess""\n"
fi

echo $paraList > parameters.txt