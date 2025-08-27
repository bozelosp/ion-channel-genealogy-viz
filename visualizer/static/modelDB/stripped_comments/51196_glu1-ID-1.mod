NEURON{
POINT_PROCESS GrC_Glu1
RANGE conccap,glu
RANGE taurec,taufacil,tauin, u0, usr 
}
UNITS{
(molar)=(1/liter)
(mM)=(millimolar)
}
PARAMETER {
taurec=800 (ms)
tauin=3 (ms)
taufacil=20(ms)
usr=0.1 
u0=0 <0,1> 
 }
ASSIGNED{
tx1(ms)
conccap 
conc01(mM)
conc02(mM)
conc03(mM)
vspr 
glu (mM)
glu1 (mM)
glu2 (mM)
glu3 (mM)
}
INITIAL {
tx1=10000000
conc01=0
conc02=0
conc03=0

conccap=0
glu=0
}
BEFORE BREAKPOINT
{ 
if (t<tx1){
glu=0
glu1=0
glu2=0
glu3=0
}
if(t>=tx1) {
glu1= conccap*(2.88(mM)+conc01)*exp((tx1-t)/0.05(ms))
glu2=conccap*(0.2(mM)+conc02)*exp((tx1-t)/0.5(ms))
glu3=conccap*(0.04(mM)+conc03)*exp((tx1-t)/1.7(ms))
glu=glu1+glu2+glu3
}
}
NET_RECEIVE (weight,Eav,R, u,tsyn (ms))
{
INITIAL
{
R=1
Eav=0
u=u0
tsyn=t}
vspr=((1-R-Eav)/taurec+(R-1)/tauin)/(1/tauin-1/taurec)
R=1+exp((tsyn-t)/taurec)*vspr+exp((tsyn-t)/tauin)*(R-1-vspr)
Eav=Eav*exp((tsyn-t)/tauin)
if (taufacil>0){
u=u*exp((tsyn-t)/taufacil)
}else {
u=usr
}
if (taufacil>0) {
u = u+usr*(1-u)
}
Eav=Eav+R*u
conccap=(u*R/usr)*weight
R=R-u*R
tsyn=t
conc01=glu1
conc02=glu2
conc03=glu3
tx1=t 
}