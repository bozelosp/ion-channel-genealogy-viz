NEURON{
POINT_PROCESS GrC_Gludif3
RANGE    glu,rPSD,h,nu,gludir,gluspill
RANGE Deff,meandist,rabs,alpha,h,Rmf  
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill 
RANGE ts1,td1,tm1
}
UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
(um)=(micron)
(nA)=(nanoamp)
PI=(pi) (1)
}
PARAMETER {
nu=2(/um2)
rabs=0(um)
h=0.02 (um)
alpha=5 
Deff=0.043 (um2/ms)
c0cleft = 8.769 (mM)
rPSD=0.11 (um) 
meandist=0.29 (um) 
Rmf=2.9(um) 
Popeak=0.678
inclugludir=1 
inclugluspill=1 
tm1=0 (ms) 
td1=0 (ms) 
ts1=0 (ms) 
 }
ASSIGNED{
   tx1(ms)
gludir (mM)
gluspill(mM)
glu (mM)
Podir
Pospill 
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
Pospill=0
}
if(t>tx1) {
UNITSOFF
gludir= (2*c0cleft*h*alpha/sqrt(4*PI*Deff*(t-tx1)))*
(1-exp(rPSD*rPSD/(4*Deff*(tx1-t))))
if(gludir>c0cleft){gludir=c0cleft}
gluspill = 2*nu*c0cleft*h*rPSD*rPSD*PI*alpha*
(1/sqrt(4*PI*Deff*(t-tx1)))*(exp(meandist*meandist/
(4*Deff*(tx1-t)))-exp(Rmf*Rmf/(4*Deff*(tx1-t))))

glu= inclugluspill*gluspill + inclugludir*gludir





Podir=(0.94*exp((tx1-t)/0.37(ms))+0.06*exp((tx1-t)/2.2(ms))
  -exp((tx1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak
Pospill=(0.39*exp((tx1-t)/2.0(ms))+0.61*exp((tx1-t)/9.1(ms))-
 exp((tx1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak
}
}
NET_RECEIVE (weight)
{
tx1=t 
}