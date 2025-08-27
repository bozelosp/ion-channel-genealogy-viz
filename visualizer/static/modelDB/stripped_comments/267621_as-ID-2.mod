UNITS
  {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (molar) = (1/liter)
  (mM) = (millimolar)
  (S) = (mho)
  (uS) = (microS)
  }
 
NEURON
  {
  SUFFIX as
  USEION k READ ek WRITE ik
  RANGE gbar, g, ik  
  }


PARAMETER
  {
  gbar = 0  (uS/mm2)
  }

ASSIGNED 
  {
  v (mV)
  ek (mV)
  minf
  taum (ms)
  hinf
  tauh (ms)
  g (uS/mm2)
  ik (mA/cm2)
  }

STATE  
       
  {
  m
  h
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m^3*h
  ik = (1e-4)*(g*(v-ek))
  }

DERIVATIVE state_change
  {
  rates(v)
  m' = (minf-m)/taum
  h' = (hinf-h)/tauh
  }

INITIAL 
  {
  rates(v)  
  m = minf
  h = hinf
  }

PROCEDURE rates(v(mV)) 
  {  
  TABLE minf, taum, hinf, tauh
  FROM -150 TO 150 WITH 301
  minf=1/(1+exp(-(v+21.0(mV))/22.8(mV)))
  taum=10.3(ms)+(5.0(ms)-10.3(ms))/(1+exp(-(v-5.6(mV))/4(mV)))
  hinf=1/(1+exp( (v+55.0(mV))/4.8(mV)))
  tauh=1/(1/(253.4(ms)*(1+exp(-(v+ 9.0(mV))/11.1(mV))))+
          1/(250.5(ms)*(1+exp( (v+92.0(mV))/16.0(mV)))))
  }