NEURON {  POINT_PROCESS GABAA }
PARAMETER {
  Cdur	= 0.3	(ms)		
  Alpha	= 10	(/ms mM)	
  Beta	= 0.16	(/ms)		
  Erev	= -80	(mV)		
}
INCLUDE "netcon.inc"