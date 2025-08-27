INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS pGpulse
  RANGE active
  RANGE g_e,g_i,g_eout,g_iout,del, dur, i0, amp_ge,amp_gi, pfreq, pdur,E_e, E_i
  USEION z READ ez WRITE iz VALENCE 1
  POINTER vcell
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
(umho) = (micromho)
}

PARAMETER {
  active = 0		
  
  del = 0	(ms)	
  dur = 10000	(ms)	
  
  E_e= 0 	  (mV)  
  E_i= -75 	  (mV)  

  amp_ge = 0	(umho)	
  amp_gi = 0	(umho)	
  pfreq = 10	(Hz)	
  pdur = 500	(ms)	
}

ASSIGNED {
  vcell	(mV)		
  iz		(nA)
  ez		(mV)
  Npulses		
  pstart	(ms)	
  pstop		(ms)	
  g_e   (umho)          
  g_i   (umho)          
  g_eout (umho)         
  g_iout (umho)      	
}

INITIAL {
  iz = 0
  g_e = 0
  g_i = 0
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

      g_e = amp_ge
      
      Npulses = Npulses + 1
      pstart = del + Npulses*1000/pfreq
      if (1000/pfreq <= pdur) { pstop = del + dur }

          } 
    else if (t >= pstop) { 
      g_e = 0
      

      pstop = pstart + pdur
    }
  }
  else { g_e = 0 
         
}

iz = g_e * (vcell - E_e) 

      g_eout = 50000.0 * g_e
      



}

PROCEDURE Update() {
  VERBATIM
    Npulses = (int) ((t-del)*pfreq/1000);
  ENDVERBATIM
  
  pstart = del + Npulses*1000/pfreq
  
  if (1000/pfreq <= pdur) { pstop = del + dur }
  else { pstop = pstart + pdur }
}