INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS pIpulse
  RANGE active
  RANGE del, dur, i0, amp, pfreq, pdur
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
}

PARAMETER {
  active = 0		
  
  del = 0	(ms)	
  dur = 10000	(ms)	
  i0 = 0	(nA)	
  amp = 0	(nA)	
  pfreq = 1	(Hz)	
  pdur = 500	(ms)	
}

ASSIGNED {
  iz		(nA)
  ez		(mV)
  Npulses		
  pstart	(ms)	
  pstop		(ms)	
}

INITIAL {
  iz = 0
  Npulses = 0
  pstart = del
  pstop = pstart + pdur

	VERBATIM
 	printf("%f \n", pstart);
	ENDVERBATIM
}

BREAKPOINT {
  if (active) {
    if ((t >= pstart) && (t < dur)) { 

	VERBATIM
 	printf("%f \n", active);
	ENDVERBATIM

      iz = -(amp + i0)
      Npulses = Npulses + 1
      pstart = del + Npulses*1000/pfreq
      if (1000/pfreq <= pdur) { pstop = del + dur }
    } 
    else if (t >= pstop) { 
      iz = -i0 
      pstop = pstart + pdur
    }
  }
  else { iz = 0 }
}

PROCEDURE Update() {
  VERBATIM
    Npulses = (int) ((t-del)*pfreq/1000);
  ENDVERBATIM
  
  pstart = del + Npulses*1000/pfreq
  
  if (1000/pfreq <= pdur) { pstop = del + dur }
  else { pstop = pstart + pdur }
}