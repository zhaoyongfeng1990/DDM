#!/bin/bash

#  compile.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

Type=$1

if [ "$Type" = "SW" ]
then
    echo "#define ISFSW" > ISFtype.h
elif [ "$Type" = "SWP" ]
then
    echo "#define ISFSWP" > ISFtype.h
elif [ "$Type" = "RTD" ]
then
    echo "#define ISFRTD" > ISFtype.h
elif [ "$Type" = "RTDP" ]
then
    echo "#define ISFRTDP" > ISFtype.h
fi

make

if [ "$Type" = "SW" ]
then
    mv ./ddm ./ddmSW
elif [ "$Type" = "SWP" ]
then
    mv ./ddm ./ddmSWP
elif [ "$Type" = "RTD" ]
then
    mv ./ddm ./ddmRTD
elif [ "$Type" = "RTDP" ]
then
    mv ./ddm ./ddmRTDP
fi
