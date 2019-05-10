#!/bin/sh

#  compileAll.sh
#  ddm
#
#  Created by Zhao Yongfeng on 15/7/20.
#  Copyright (c) 2015年 ZYF. All rights reserved.

echo "#define ISFRTDP" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDP_multisets

echo "#define ISFRTDPfix" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDPfix_multisets

echo "#define ISFRTDPTT" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDPTT_multisets

echo "#define ISFRTDPTTfix" > ISFtype.h
make
mv ./ddm ./bin/ddmRTDPTTfix_multisets