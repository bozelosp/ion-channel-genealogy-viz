INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS cIpulse
  RANGE active
  RANGE del, dur, amp 
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
}

PARAMETER {
  active = 0		
  
  del = 0	(ms)	
  dur = 1e+09	(ms)	
  amp = 0	(nA)	
}

ASSIGNED {
  iz	(nA)
  ez	(mV)
}

INITIAL {
  iz = 0
}

BREAKPOINT {
  if (active)  {
    if ((t >= del) && (t < dur)) { iz = -amp }
    else { iz = 0 }
  }
  else { iz = 0 }
}

PROCEDURE Update() {
  
}