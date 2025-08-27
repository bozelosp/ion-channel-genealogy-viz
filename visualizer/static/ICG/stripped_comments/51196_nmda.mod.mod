NEURON{
POINT_PROCESS  GrC_NMDA
POINTER pNp, pNd, pglu
RANGE   inmda, icanmda
RANGE   gnmda, gbarnmda,freqdel 
RANGE   Erev 
NONSPECIFIC_CURRENT  inmda 
}
UNITS{
(nA)=(nanoamp)
(mA)=(milliamp)
(mV)=(millivolt)
(umho)=(micromho)
(molar)=(1/liter)
(mM)=(millimolar)
}
PARAMETER {
Erev = 0 (mV)
gbarnmda = 9.48e-3 (umho)
alpha_vspom=-0.062 (/mV)

v0_block=10 (mV)
cnm_4_3 = 10 (/mM-ms)
cnm_3_2 = 5 (/mM-ms)
cnm_1_2 = 0.2 (/ms)
cnm_2_1 = 0.06 (/ms)
cnm_2_3 = 0.0134 (/ms)
cnm_2_5 = 0.01 (/ms)
cnm_3_4 = 0.0067 (/ms)
cnm_5_2= 0.008 (/ms)
freqdel=  0.579 (/nA)
}
ASSIGNED{
v (mV)
inmda (nA)
icanmda (nA)
gnmda (umho)
pNp (nA) 
pNd (nA) 
pglu (mM) 
celsius(degC)
}
STATE{
St1nmda
St2nmda
St3nmda
St4nmda
St5nmda
  }
INITIAL {
St1nmda=0
St2nmda=0
St3nmda=0
St4nmda=1
St5nmda=0
}

BREAKPOINT
{
SOLVE states METHOD derivimplicit
gnmda = gbarnmda*St1nmda *(1+freqdel*(pNp-pNd))
inmda =gnmda*vspom(v)*(v-Erev)
icanmda=0.1*gnmda*(v-26)*vspom(v)
}

DERIVATIVE states {
LOCAL Q10
Q10 = 2^((celsius-20(degC))/10(degC))
St1nmda' = Q10*( - cnm_1_2*St1nmda + cnm_2_1*St2nmda)
St2nmda' = Q10*(-(cnm_2_1 + cnm_2_3  + cnm_2_5) * St2nmda + 
cnm_1_2*St1nmda + cnm_3_2*pglu*St3nmda + cnm_5_2*St5nmda)
St3nmda' = Q10*(-(cnm_3_2*pglu+cnm_3_4) * St3nmda +
cnm_4_3 *pglu* St4nmda + cnm_2_3* St2nmda)
St4nmda' = Q10*(-cnm_4_3*pglu*St4nmda+cnm_3_4*St3nmda)
St5nmda' = Q10*(-cnm_5_2*St5nmda+ cnm_2_5*St2nmda)
}

FUNCTION vspom (v(mV))( ){
vspom=1./(1.+0.2801*1.2*exp(alpha_vspom*(v-v0_block)))
}