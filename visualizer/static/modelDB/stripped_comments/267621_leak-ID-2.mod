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
  SUFFIX leak
  NONSPECIFIC_CURRENT i
  RANGE gbar, g, i, e  
  }


PARAMETER
  {
  gbar = 0 (uS/mm2)
  e = -50 (mV)
  }

ASSIGNED
  {
  v (mV)
  g (uS/mm2)
  i (mA/cm2)
  }

BREAKPOINT
  {
  g=gbar
  i = (1e-4)*(g*(v-e))
  }