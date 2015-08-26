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
    if [ "$Type" = "RD" ]
    then
        echo "#define ISFRD" > ISFtype.h
    elif [ "$Type" = "RDP" ]
    then
        echo "#define ISFRDP" > ISFtype.h
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

if [ "$Type" = "RD" ]
then
    mv ./ddm ./bin/ddmRD
elif [ "$Type" = "RDP" ]
then
    mv ./ddm ./bin/ddmRDP
elif [ "$Type" = "RTD" ]
then
    mv ./ddm ./bin/ddmRTD
elif [ "$Type" = "RTDP" ]
then
    mv ./ddm ./bin/ddmRTDP
elif [ "$Type" = "RTDPfix" ]
then
    mv ./ddm ./bin/ddmRTDPfix
fi

