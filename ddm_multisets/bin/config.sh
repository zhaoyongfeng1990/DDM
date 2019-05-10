#!/bin/sh

#  config.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

omp_num=1
nqcurve=2

qincre=30

dimx=512
dimy=512
numOfSeq=4500
numOfDiff=4000
dx=2.6

qmin=0.01
qstep=0.01
dt=0.01
timeWindow=10

maxIter=200
alphaGuess=0.8
DGuess=0.2
vbarGuess=13
lambdaGuess=0.5
ZGuess=15
sigmaGuess=3
TTGuess=0.05

if [ $nqcurve = 1 ]
then
    paraList="$omp_num""\n""$nqcurve""\n""$dimx""\n""$dimy""\n""$numOfSeq""\n""$numOfDiff""\n""$dx""\n""$qmin""\n""$qstep""\n""$dt""\n""$timeWindow""\n""$maxIter""\n""$alphaGuess""\n""$DGuess""\n""$vbarGuess""\n""$lambdaGuess""\n""$ZGuess""\n""$sigmaGuess""\n""$TTGuess""\n"
else
    paraList="$omp_num""\n""$nqcurve""\n""$qincre""\n""$dimx""\n""$dimy""\n""$numOfSeq""\n""$numOfDiff""\n""$dx""\n""$qmin""\n""$qstep""\n""$dt""\n""$timeWindow""\n""$maxIter""\n""$alphaGuess""\n""$DGuess""\n""$vbarGuess""\n""$lambdaGuess""\n""$ZGuess""\n""$sigmaGuess""\n""$TTGuess""\n"
fi

echo $paraList > parameters.txt
