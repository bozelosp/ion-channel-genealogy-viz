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
  SUFFIX kd
  USEION k READ ek WRITE ik
  RANGE gbar, g, ik
  }


PARAMETER
  {
  gbar = 0 (uS/mm2)
  }

ASSIGNED
  {
  v (mV)
  ek (mV)
  minf
  taum (ms)
  g (uS/mm2)
  ik (mA/cm2)
  }

STATE
  {
  m
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m^4
  ik = (1e-4)*(g*(v-ek))
  }

DERIVATIVE state_change
  {
  rates(v)  
  m' = (minf-m)/taum
  }

INITIAL 
  {
  rates(v)  
  m = minf
  }

PROCEDURE rates(v(mV)) 
  {  
  TABLE minf, taum
  FROM -150 TO 150 WITH 301
  minf=1/(1+exp(-(v+25(mV))/17(mV)))
  taum=120(ms)+(6.2(ms)-120(ms))/(1+exp(-(v+46.1(mV))/18.1(mV)))
  }