NEURON {  POINT_PROCESS GABAb }

PARAMETER {
  Cdur	= 85	(ms)		
  Alpha	= 0.016	(/ms mM)	
  Beta	= 0.0047 (/ms)		
  Erev	= -90	(mV)		
  DELAY = 0
  Deadtime = 1. (ms)		
  GMAX	= 1.0	(uS)		
  Thresh = -70
}

INCLUDE "sns.inc"

    
    BREAKPOINT {
      if (v>Thresh) {
        g = g*((v/Thresh)*(v/Thresh)*(v/Thresh))
        i = g*(v - Erev)
      }
    }