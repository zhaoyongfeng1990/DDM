#!/bin/bash

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
    if [ "$Type" = "RTDP" ]
    then
        echo "#define ISFRTDP" > ISFtype.h
    elif [ "$Type" = "RTDPfix" ]
    then
        echo "#define ISFRTDPfix" > ISFtype.h
    elif [ "$Type" = "RTDPTT" ]
    then
        echo "#define ISFRTDPTT" > ISFtype.h
    elif [ "$Type" = "RTDPTTfix" ]
    then
        echo "#define ISFRTDPTTfix" > ISFtype.h
    fi
fi

make

if [ "$Type" = "RTDP" ]
then
    mv ./ddm ./bin/ddmRTDP_multisets
elif [ "$Type" = "RTDPfix" ]
then
    mv ./ddm ./bin/ddmRTDPfix_multisets
elif [ "$Type" = "RTDPTT" ]
then
    mv ./ddm ./bin/ddmRTDPTT_multisets
elif [ "$Type" = "RTDPTTfix" ]
then
    mv ./ddm ./bin/ddmRTDPTTfix_multisets
fi

