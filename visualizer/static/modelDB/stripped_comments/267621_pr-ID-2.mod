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
  SUFFIX pr
  NONSPECIFIC_CURRENT i
  RANGE e, gbar, i,g, theta 
  }


PARAMETER
  {
  gbar = 0 (uS/mm2)
  e = -10 (mV)
  theta = -45 (mV)
  }

ASSIGNED
  {
  v (mV)
  minf
  g (uS/mm2)
  i (mA/cm2)
  }

STATE
  {
  m
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m
  i = (1e-4)*(g*(v-e)) 
  }

DERIVATIVE state_change
  {
  rates(v)
  m' = (minf-m)/(6(ms))
  }

INITIAL 
  {
  rates(v)    
  m = minf
  }

PROCEDURE rates(v(mV)) 
  {


  minf=1/(1+exp(-(v-theta)/5(mV)))
  }