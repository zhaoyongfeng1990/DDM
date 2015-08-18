#!/bin/sh

#  compile.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

Type=$1

ISF=$(cat ./ISFtype.h)

OldType=$(echo ${ISF:11})

if [ "$Type" != "$OldType" ]
then
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
    elif [ "$Type" = "RTDPfix" ]
    then
        echo "#define ISFRTDPfix" > ISFtype.h
    fi
fi

make

if [ "$Type" = "SW" ]
then
    mv ./ddm ./bin/ddmSW
elif [ "$Type" = "SWP" ]
then
    mv ./ddm ./bin/ddmSWP
elif [ "$Type" = "RTD" ]
then
    mv ./ddm ./bin/ddmRTD
elif [ "$Type" = "RTDP" ]
then
    mv ./ddm ./bin/ddmRTDP
elif [ "$Type" = "RTDPfix" ]
then
    mv ./ddm ./bin/ddmRTDP
fi

