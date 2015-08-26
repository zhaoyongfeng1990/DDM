#!/bin/sh

#  compileAll.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

echo "#define ISFRD" > ISFtype.h
make
mv ./ddm ./bin/ddmRD

echo "#define ISFRDP" > ISFtype.h
make
mv ./ddm ./bin/ddmRDP

echo "#define ISFRTD" > ISFtype.h
make
mv ./ddm ./bin/ddmRTD

echo "#define ISFRTDP" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDP

echo "#define ISFRTDPfix" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDPfix
