#!/bin/sh

#  compileAll.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

echo "#define ISFSW" > ISFtype.h
make
mv ./ddm ./bin/ddmSW

echo "#define ISFSWP" > ISFtype.h
make
mv ./ddm ./bin/ddmSWP

echo "#define ISFRTD" > ISFtype.h
make
mv ./ddm ./bin/ddmRTD

echo "#define ISFRTDP" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDP

echo "#define ISFRTDPfix" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDPfix
