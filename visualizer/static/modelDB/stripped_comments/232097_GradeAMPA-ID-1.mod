NEURON {
  POINT_PROCESS gradAMPA
  RANGE gmax, g, i, alpha, beta, thetasyn,e, sigma
  
  NONSPECIFIC_CURRENT i
  POINTER vpre
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  gmax  = 1e-3 (uS)
  alpha = 1 (/ms)      
  beta  = 0.1818 (/ms) 
  e =  0	  (mV)     
  thetasyn =  0 (mV)   
  sigma    = 2  
}

ASSIGNED { vpre (mV) v (mV) i (nA)  g (uS)}

STATE { s }

INITIAL {
  s =  alpha*F(vpre)/(alpha*F(vpre)+beta)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  g = gmax * s
  i = g*(v - e)
}

DERIVATIVE state {
  s' = alpha*F(vpre)*(1-s) - beta*s
}

FUNCTION F (v1 (mV)) {
  F = 1/(1 + exp(-(v1-thetasyn)/sigma))
}