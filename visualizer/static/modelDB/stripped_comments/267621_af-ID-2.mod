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
  SUFFIX af
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
  m' = (minf-m)/(3(ms))
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
  TABLE minf, hinf, tauh
  FROM -150 TO 150 WITH 301
  minf=1/(1+exp(-(v+14.5(mV))/18.1(mV)))
  hinf=1/(1+exp( (v+68.1(mV))/4.5(mV)))
  tauh=119.4(ms)+(19.3(ms)-119.4(ms))/(1+exp(-(v+1.8(mV))/4.0(mV)))
  }