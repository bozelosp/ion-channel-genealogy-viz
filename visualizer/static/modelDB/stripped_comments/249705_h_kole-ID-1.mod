UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
}

PARAMETER {
  v (mV)
  erev=-45  		(mV) 	
  gbar=0.00015 	(S/cm2)	
  q10 = 2.2
  ascale = 0.00643
  bscale = 0.193
  ashift = 154.9
  aslope = 11.9
  bslope = 33.1
}

NEURON {
 THREADSAFE
  SUFFIX ih
  NONSPECIFIC_CURRENT i
  RANGE i,gbar,ascale,bscale,ashift,aslope,bslope
}

STATE {
  m
}

ASSIGNED {
  i (mA/cm2)
}

INITIAL { LOCAL a,b
  a = alpha(v)
  b = beta(v)
  m = a / (a + b)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  i = gbar*m*(v-erev)
}


FUNCTION alpha(v(mV)) {
  alpha = ascale*(v+ashift)/(exp((v+ashift)/aslope)-1)  
  
}

FUNCTION beta(v(mV)) {
  beta = bscale*exp(v/bslope)
}

DERIVATIVE state {
  m' = (1-m)*alpha(v) - m*beta(v)
}