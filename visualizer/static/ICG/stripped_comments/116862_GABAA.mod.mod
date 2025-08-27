NEURON {  POINT_PROCESS GABAA }
PARAMETER {
  Cdur	= 1.0	(ms)		
  Alpha	= 0.53	(/ms mM)	
  Beta	= 0.18	(/ms)		
  Erev	= -80	(mV)		
}
INCLUDE "netcon.inc"