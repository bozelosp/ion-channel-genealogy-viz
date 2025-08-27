NEURON{
POINT_PROCESS GrC_Gludif2
RANGE    PRE, glu,rPSD,nu,gludir,gluspill
RANGE Deff,meandist,rabs,Rmf 
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill, Podir1
RANGE td1,tm1,ts1
}

UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
(um)=(micron)
(nA)=(nanoamp)
PI=(pi)  (1)
}

PARAMETER {
alpha=1   
nu=0.33(/um2)  
rabs=0 (um)     
Deff=0.2 (um2/ms)
c0cleft = 8.769 (mM) 
rPSD=0.11 (um)   
meandist=0.23 (um)
Rmf=2.8 (um) 
Popeak=0.634 
inclugludir=1
inclugluspill=1
td1=0 (ms) 
tm1=0 (ms) 
ts1=0 (ms) 
}
ASSIGNED{
Podir
Podir1
Pospill 
tx1(ms)
gludir (mM)
gluspill(mM)
vspr
glu (mM)
}
INITIAL {
tx1=10000000
glu=0 
gludir=0
gluspill=0
}
BREAKPOINT
{
at_time(tx1)
if (t<=tx1){
glu=0
gludir=0
gluspill=0
Podir=0
Podir1=0
Pospill=0
}
if(t>tx1) {
gludir= c0cleft*(1-exp(rPSD*rPSD/(4*Deff*(tx1-t))))
if(gludir>c0cleft){gludir=c0cleft}
gluspill=PI*nu*c0cleft*rPSD*rPSD*
(exp(meandist*meandist/(4*Deff*(tx1-t)))-exp(Rmf*Rmf/(4*Deff*(tx1-t))))
 glu= inclugludir*gludir  +inclugluspill*gluspill




Podir=(0.94*exp((tx1-t)/0.37(ms))+0.06*exp((tx1-t)/2.2(ms))
  -exp((tx1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak
Podir1=(0.94*exp((tx1-t)/0.3(ms))+0.06*exp((tx1-t)/3.1(ms))
  -exp((tx1-t)/0.12(ms)))/0.35*Popeak 
Pospill=(0.39*exp((tx1-t)/2.0(ms))+0.61*exp((tx1-t)/9.1(ms))-
 exp((tx1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak
}
}
NET_RECEIVE (weight)
{
tx1=t 
}