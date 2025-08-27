NEURON {POINT_PROCESS AMPA}

PARAMETER {
  Cdur	= 0.3	(ms)		
  Alpha	= 0.47	(/ms mM)	
  Beta	= 0.18	(/ms)		
  Erev	= 0	(mV)		
}
INCLUDE "netcon.inc"