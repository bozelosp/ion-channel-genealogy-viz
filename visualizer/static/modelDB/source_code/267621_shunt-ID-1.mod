: A shunt current

NEURON
  {
  POINT_PROCESS Shunt
  NONSPECIFIC_CURRENT I
  RANGE I, e, G
  }

PARAMETER
  {
  G = 0 (micromho)
  e = 0 (millivolt) 
  }

ASSIGNED
  {
  I  (nanoamp)
  v  (millivolt)  
  }

BREAKPOINT
  { 
  I = G * (v - e)
  }

