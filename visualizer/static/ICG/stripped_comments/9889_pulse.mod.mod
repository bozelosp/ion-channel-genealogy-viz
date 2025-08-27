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
  DELAY = 0                     
  dur	= 0	(ms)            
  amp   = 0     (nA)            
}

ASSIGNED {
  Aon
  i 		(nA)
  dt
}

INITIAL {
  i = 0
}

BREAKPOINT {
  i = Aon
}


NET_RECEIVE(weight, on, nspike) {
  if (flag == 0 && !on) { 
    nspike = nspike + 1 
    Aon = -amp
    
    net_send(dur, nspike)
  }
  if (flag == nspike) { 
    Aon = 0
    on = 0
  }
}