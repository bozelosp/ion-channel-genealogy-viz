NEURON {  POINT_PROCESS GABAa }

PARAMETER {
  Cdur	= 1.08	(ms)		
  Alpha	= 1.	(/ms mM)	
  Beta	= 0.5	(/ms)		
  Erev	= -75	(mV)		
  DELAY = 0
  Deadtime = 0 (ms)		
  GMAX	= 1.0	(uS)		
}

INCLUDE "sns.inc"