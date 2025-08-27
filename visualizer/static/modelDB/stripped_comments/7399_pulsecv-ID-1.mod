INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PULSE
  NONSPECIFIC_CURRENT i
  RANGE amp, dur
}
 
UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  dur	= 0	(ms)            
  amp   = 0     (nA)            
}

ASSIGNED {
  on
  i 		(nA)
  dt
}

INITIAL {
  on = 0
  i = 0
}

BREAKPOINT {
  if (on==1) { 
    i = -amp
  } else {
    i = 0
  }
}

NET_RECEIVE(weight) {   
  if (flag == 0) { 
    on = 1
    net_send(dur,-1) 
  }
  if (flag == -1) { 
    on = 0
  }
}