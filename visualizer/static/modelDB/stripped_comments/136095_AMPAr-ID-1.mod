NEURON {  POINT_PROCESS AMPAr }

PARAMETER {
  Cdur	= 1	(ms)		
  Alpha	= 1.	(/ms mM)	
  Beta	= 0.5	(/ms)		
  Erev	= 0	(mV)		
}

INCLUDE "netrand.inc"