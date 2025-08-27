NEURON {  POINT_PROCESS GABAA }









PARAMETER {
  Cdur	= 1.0	(ms)		
  Alpha	= 0.53	(/ms mM)	
  Beta	= 0.18	(/ms)		
  Erev	= -75	(mV)		
  Deadtime = 1.	(ms)		
  GMAX = 1	(uS)		
  DELAY = 2
}
INCLUDE "netcon.inc"