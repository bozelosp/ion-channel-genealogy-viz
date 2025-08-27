NEURON {  POINT_PROCESS GABAar }

PARAMETER {
  Cdur	= 1.08	(ms)		
  Alpha	= 1.	(/ms mM)	
  Beta	= 0.5	(/ms)		
  Erev	= -70	(mV)		
}

INCLUDE "netrand.inc"