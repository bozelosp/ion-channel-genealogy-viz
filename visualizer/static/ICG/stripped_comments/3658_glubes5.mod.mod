NEURON{
POINT_PROCESS GrC_Glubes5
RANGE  glu,rPSD,rabs,nu, gluspill,gludir,Rmf
RANGE Deff,meandist,rabs 
RANGE inclugludir,inclugluspill, Popeak,alpha,Podir,Pospill
RANGE tm1,td1,ts1
}
UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
(um)=(micron)
(nA)=(nanoamp)
}
CONSTANT {
PI=3.1415927
}
PARAMETER { 
alpha=1 
Deff=0.2 (um2/ms) a
nu=0.33(1/um2)  
rabs= 4.1 (um) 
c0cleft = 8.769 (mM)
rPSD=0.11 (um)
meandist=0.23 (um)
Rmf=2.8(um) 
Popeak=0.634 
inclugludir=1 
inclugluspill=1 
tm1=0 (ms) 
td1=0 (ms) 
ts1=0 (ms) 
}
VERBATIM
int i;
double l[100000];
extern float bessj1();
extern float bessj0();
ENDVERBATIM
ASSIGNED{
Podir
Pospill 
tx1(ms)
gludir (mM)
gluspill(mM)
vspr
glu (mM)
sum (um)
sum0 (um)
sum1(um2)
sum01(um2)
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
Podir=0
Pospill=0
gludir=0
gluspill=0
}
if(t>tx1) {
VERBATIM
l[0]=0;
l[1]=3.8317;l[2]=7.01558667;l[3]=10.1737;
sum=0; i=1; 
do 
{if (i>=4) 
l[i]=PI*(4*i+1)/4;
sum0=sum;
sum =sum+bessj1((l[i]/rabs)*rPSD)/((l[i]/rabs)*bessj0(l[i])*bessj0(l[i]))
* exp((l[i]/rabs)*(l[i]/rabs)*Deff*(tx1-t));
 i++; }
while (fabs(sum-sum0)>=0.001);
ENDVERBATIM
gludir= c0cleft*2*rPSD*((rPSD/2)+sum)/(rabs*rabs)
if(gludir>c0cleft){gludir=c0cleft}
VERBATIM
sum1=0; i=1;
do
{if (i>=4) 
l[i]=PI*(4*i+1)/4;
sum01=sum1;
sum1=sum1+(Rmf*bessj1((l[i]/rabs)*Rmf)-meandist 
* bessj1((l[i]/rabs) *meandist))/
((l[i]/rabs)*bessj0(l[i])*bessj0(l[i]))*exp((l[i]/rabs)*(l[i]/rabs)
*Deff*(tx1-t));
i++;}
while (fabs(sum1-sum01)>=0.0001);
ENDVERBATIM
gluspill=2*PI*nu*c0cleft*rPSD*rPSD*((Rmf*Rmf-meandist*meandist)/2+sum1)/(rabs*rabs)
glu= inclugludir*gludir  +inclugluspill*gluspill


Podir=(0.94*exp((tx1-td1-t)/0.37(ms))+0.06*exp((tx1-td1-t)/2.2(ms))
  -exp((tx1-td1-t)/0.199(ms)))/0.249*(0.43/0.484)*Popeak
Pospill=(0.39*exp((tx1-ts1-t)/2.0(ms))+0.61*exp((tx1-ts1-t)/9.1(ms))-
 exp((tx1-ts1-t)/0.44(ms)))/0.682*(0.125/0.484)*Popeak 
}
}
NET_RECEIVE (weight)
{
tx1=t 
}