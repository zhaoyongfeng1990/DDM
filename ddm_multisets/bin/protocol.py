# -*- coding: utf-8 -*-
import numpy as np
import sys
import subprocess as subp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
	print('Not enough arguments')
	sys.exit()

omp_num=16
nqcurve=2

qincre=30

dimx=512
dimy=512
numOfSeq=4500
numOfDiff=4000
dx=1.625

qmin=0.01510381083456
qstep=0.00755190541728
dt=0.01
timeWindow=10

maxIter1=1000
maxIter2=200
alphaGuess=0.794639812500000
DGuess=0.418209969696970
vbarGuess=12.9807261538462
lambdaGuess=0.5
ZGuess=3
sigmaGuess=3.25101489339998
TTGuess=0.05

#xrange=[0,q[lq-1]]
xrange=[0,2]
alpharange=[0,1]
vrange=[0,15]
sigmarange=[0,10]
Drange=[0,1]
lambdarange=[0,10]
Zrange=[0,20]

exactAlpha=0.8
exactV=13
exactSigma=3.25
exactZ=15
exactD=0.4
exactLambda=1

paraRDP=[omp_num, 1, dimx, dimy, numOfSeq, numOfDiff, dx, qmin, qstep, dt, timeWindow, maxIter1, alphaGuess, DGuess, vbarGuess, lambdaGuess, ZGuess, sigmaGuess, TTGuess]
parafile=open('parameters.txt','w')
for i in paraRDP:
	parafile.write(str(i))
	parafile.write('\n')
parafile.close()

#os.system("./ddmRDP "+sys.argv[1])
#os.system("mkdir RDPfit")
#subp.check_call("./ddmRDP "+sys.argv[1], shell=True)
subp.check_call("mkdir RDPfit", shell=True)



# Begin first RDP plot
q=np.loadtxt('q.txt')
fitpara=np.loadtxt('fitparafile.txt')
status=np.loadtxt('statusq.txt')
q=q[status==0]
fitpara=fitpara[status==0,:]
lq=q.size

f=plt.figure()
f.figsize=(20,15)
bigFig=f.add_subplot(111)
alphaFig=f.add_subplot(231)
vFig=f.add_subplot(232)
sigmaFig=f.add_subplot(233)
ZFig=f.add_subplot(234)
DFig=f.add_subplot(235)

bigFig.spines['top'].set_color('none')
bigFig.spines['bottom'].set_color('none')
bigFig.spines['left'].set_color('none')
bigFig.spines['right'].set_color('none')
bigFig.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

alphaFig.axhline(y=exactAlpha,xmin=0,xmax=xrange[1],color='r')
alphaFig.plot(q, fitpara[:,0],'.', markersize=2)
alphaFig.set_xlim(xrange)
alphaFig.set_ylim(alpharange)
alphaFig.set_ylabel(r'$\alpha$')

vFig.axhline(y=exactV,xmin=0,xmax=xrange[1],color='r')
vFig.plot(q, fitpara[:,2],'.', markersize=2)
vFig.set_xlim(xrange)
vFig.set_ylim(vrange)
vFig.set_ylabel(r'$v$ ($\mu\mathrm{m}/\mathrm{s}$)')

DFig.axhline(y=exactD,xmin=0,xmax=xrange[1],color='r')
DFig.plot(q, fitpara[:,1],'.', markersize=2)
DFig.set_xlim(xrange)
DFig.set_ylim(Drange)
DFig.set_ylabel(r'$D$ ($\mu\mathrm{m}^2/\mathrm{s}$)')

sigmaFig.axhline(y=exactSigma,xmin=0,xmax=xrange[1],color='r')
sigmaFig.plot(q, fitpara[:,2]/np.sqrt(fitpara[:,3]+1),'.', markersize=2)
sigmaFig.set_xlim(xrange)
sigmaFig.set_ylim(vrange)
sigmaFig.set_ylabel(r'$\sigma$ ($\mu\mathrm{m}/\mathrm{s}$)')

ZFig.axhline(y=exactZ,xmin=0,xmax=xrange[1],color='r')
ZFig.plot(q, fitpara[:,3],'.', markersize=2)
ZFig.set_xlim(xrange)
ZFig.set_ylim(Zrange)
ZFig.set_ylabel(r'$Z$')

bigFig.set_xlabel(r'$q/\mu\mathrm{m}^{-1}$')
bigFig.set_title('RDP fitting result of '+sys.argv[1])

plt.tight_layout(pad=0.2, w_pad=0.1, h_pad=0.2)

f.savefig('./RDPfit/RDPfit.png',dpi=300)

#os.system('cp ./*.txt ./RDPfit/')
subp.check_call('cp ./*.txt ./RDPfit/', shell=True)

# extract estimation
avePara=[]
for i in range(0,4):
    plist=fitpara[(q>1) & (q<1.5) ,i]
    while True:   
        meanP=np.mean(plist)
        stdP=np.std(plist)
        lp=plist.size
        plist=plist[np.abs(plist-meanP)<3*stdP]
        if plist.size==lp:
            avePara.append(meanP)
            break

sigmaGuess=avePara[2]/np.sqrt(avePara[3]+1)

paraRTDP=[omp_num, nqcurve, qincre, dimx, dimy, numOfSeq, numOfDiff, dx, qmin, qstep, dt, timeWindow, maxIter2, avePara[0], avePara[1], avePara[2], lambdaGuess, ZGuess, sigmaGuess, TTGuess]
parafile=open('parameters.txt','w')
for i in paraRTDP:
	parafile.write(str(i))
	parafile.write('\n')
parafile.close()

#os.system("./ddmRTDP recover")
#os.system("mkdir RTDPfit")
subp.check_call("./ddmRTDP recover", shell=True)
subp.check_call("mkdir RTDPfit", shell=True)

# Begin first RTDP plot
q=np.loadtxt('q.txt')
fitpara=np.loadtxt('fitparafile.txt')
status=np.loadtxt('statusq.txt')
q=q[status==0]
fitpara=fitpara[status==0,:]
lq=q.size

f=plt.figure()
f.figsize=(20,15)
bigFig=f.add_subplot(111)
alphaFig=f.add_subplot(231)
vFig=f.add_subplot(232)
sigmaFig=f.add_subplot(233)
ZFig=f.add_subplot(234)
lambdaFig=f.add_subplot(235)
DFig=f.add_subplot(236)

bigFig.spines['top'].set_color('none')
bigFig.spines['bottom'].set_color('none')
bigFig.spines['left'].set_color('none')
bigFig.spines['right'].set_color('none')
bigFig.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

alphaFig.axhline(y=exactAlpha,xmin=0,xmax=xrange[1],color='r')
alphaFig.plot(q, fitpara[:,0],'.', markersize=2)
alphaFig.set_xlim(xrange)
alphaFig.set_ylim(alpharange)
alphaFig.set_ylabel(r'$\alpha$')

vFig.axhline(y=exactV,xmin=0,xmax=xrange[1],color='r')
vFig.plot(q, fitpara[:,1],'.', markersize=2)
vFig.set_xlim(xrange)
vFig.set_ylim(vrange)
vFig.set_ylabel(r'$v$ ($\mu\mathrm{m}/\mathrm{s}$)')

DFig.axhline(y=exactD,xmin=0,xmax=xrange[1],color='r')
DFig.plot(q, fitpara[:,4],'.', markersize=2)
DFig.set_xlim(xrange)
DFig.set_ylim(Drange)
DFig.set_ylabel(r'$D$ ($\mu\mathrm{m}^2/\mathrm{s}$)')

sigmaFig.axhline(y=exactSigma,xmin=0,xmax=xrange[1],color='r')
sigmaFig.plot(q, fitpara[:,2],'.', markersize=2)
sigmaFig.set_xlim(xrange)
sigmaFig.set_ylim(vrange)
sigmaFig.set_ylabel(r'$\sigma$ ($\mu\mathrm{m}/\mathrm{s}$)')

ZFig.axhline(y=exactZ,xmin=0,xmax=xrange[1],color='r')
ZFig.plot(q, fitpara[:,1]**2/fitpara[:,2]**2-1,'.', markersize=2)
ZFig.set_xlim(xrange)
ZFig.set_ylim(Zrange)
ZFig.set_ylabel(r'$Z$')

lambdaFig.axhline(y=exactLambda,xmin=0,xmax=xrange[1],color='r')
#lambdaFig.axhline(y=5,xmin=0,xmax=xrange[1],color='r')
#lambdaFig.axhline(y=3,xmin=0,xmax=xrange[1],color='g')
lambdaFig.plot(q, fitpara[:,3],'.', markersize=2)
lambdaFig.set_xlim(xrange)
lambdaFig.set_ylim(lambdarange)
lambdaFig.set_ylabel(r'$\lambda$ ($\mathrm{s}^{-1}$)')

bigFig.set_xlabel(r'$q/\mu\mathrm{m}^{-1}$')
bigFig.set_title('RTDP fitting result of '+sys.argv[1])

plt.tight_layout(pad=0.2, w_pad=0.1, h_pad=0.2)

f.savefig('./RTDPfit/RTDPfit.png',dpi=300)
#os.system('cp ./*.txt ./RTDPfit/')
subp.check_call('cp ./*.txt ./RTDPfit/', shell=True)

# extract estimation
avePara=[]
for i in range(0,5):
    plist=fitpara[(q>1) & (q<1.5) ,i]
    while True:   
        meanP=np.mean(plist)
        stdP=np.std(plist)
        lp=plist.size
        plist=plist[np.abs(plist-meanP)<3*stdP]
        if plist.size==lp:
            avePara.append(meanP)
            break

paraRTDPfix=[omp_num, nqcurve, qincre, dimx, dimy, numOfSeq, numOfDiff, dx, qmin, qstep, dt, timeWindow, maxIter2, avePara[0], avePara[4], avePara[1], avePara[3], ZGuess, avePara[2], TTGuess]
parafile=open('parameters.txt','w')
for i in paraRTDPfix:
	parafile.write(str(i))
	parafile.write('\n')
parafile.close()

#os.system("./ddmRTDPfix recover")
#os.system("mkdir lambdaFit")
subp.check_call("./ddmRTDPfix recover", shell=True)
subp.check_call("mkdir lambdaFit", shell=True)

# Begin RTDPfix plot
q=np.loadtxt('q.txt')
fitpara=np.loadtxt('fitparafile.txt')
status=np.loadtxt('statusq.txt')
q=q[status==0]
fitpara=fitpara[status==0,:]
lq=q.size

plt.figure()
plt.figsize=(20,15)

plt.axhline(y=exactLambda,xmin=0,xmax=xrange[1],color='r')
#plt.axhline(y=5,xmin=0,xmax=xrange[1],color='r')
#plt.axhline(y=3,xmin=0,xmax=xrange[1],color='g')
plt.plot(q, fitpara[:,0],'.', markersize=2)
plt.xlim(xrange)
plt.ylim(lambdarange)
plt.ylabel(r'$\lambda$ ($\mathrm{s}^{-1}$)')
plt.xlabel(r'$q/\mu\mathrm{m}^{-1}$')
plt.title(r'$\lambda$ fitting result of '+sys.argv[1])

plt.tight_layout(pad=0.2, w_pad=0.1, h_pad=0.2)

plt.savefig('./lambdaFit/lambdaFit.png',dpi=300)
#os.system('cp ./*.txt ./lambdaFit/')
subp.check_call('cp ./*.txt ./lambdaFit/', shell=True)