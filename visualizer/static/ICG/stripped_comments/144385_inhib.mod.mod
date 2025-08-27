UNITS {
  (mV) = (millivolt)
  (mA) = (milliamp)
  (S) = (siemens)
  (nA) = (nanoamp)
  (uS) = (microsiemens)
}
NEURON {
  SUFFIX inhib
  POINT_PROCESS inhib
  NONSPECIFIC_CURRENT i
  RANGE g, e
}

PARAMETER {
  g = 0.41e-3 (uS) <0,1e9> 
  e = -70 (mV) 
}

ASSIGNED {v (mV) i (nA)}
BREAKPOINT {
  i = g*(v - e)
}