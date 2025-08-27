NEURON {  POINT_PROCESS GABAb }

PARAMETER {
  Cdur	= 85	(ms)		
  Alpha	= 0.016	(/ms mM)	
  Beta	= 0.0047 (/ms)		
  Erev	= -90	(mV)		
}

INCLUDE "netcon.inc"