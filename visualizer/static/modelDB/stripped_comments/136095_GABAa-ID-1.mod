NEURON {  POINT_PROCESS GABAa }

PARAMETER {
  Cdur	= 1.08	(ms)		
  Alpha	= 1.	(/ms mM)	
  Beta	= 0.5	(/ms)		
  Erev	= -70	(mV)		
}

INCLUDE "netcon.inc"