INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PULSE
  NONSPECIFIC_CURRENT i
  RANGE amp, dur
}
 
INCLUDE "snsarr.inc"  

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  dur	= 0	(ms)            
  amp   = 0     (nA)            
  Cdur  = 0                     
  Deadtime = 0	(ms)		
  GMAX = 1	(umho)		
  DELAY = 0     (ms)
}

ASSIGNED {
  i 		(nA)
  dt
}

INITIAL {
  i = 0
  if (nsyn > 0) {
    initq()   
  } 
}

BREAKPOINT {
  if (nsyn>0) { 
  VERBATIM 
  static int ii,who;
  static QueU *pqueu;
  static SynS *ppst;

  pqueu = (QueU *)((unsigned long) queu);

  while (t >= pqueu[(int)begsyn].time) { 
    i = -amp;
    popqh1(dur);		
  }

  while (t >= pqueu[(int)endsyn].time) { 
    i = 0;
    popqh2();  
  }
  ENDVERBATIM
}
}