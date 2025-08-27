NEURON {  POINT_PROCESS GABALOW }









PARAMETER {
  Cdur	= 1.0	(ms)		
  Alpha	= 0.53	(/ms mM)	
  Beta	= 0.18	(/ms)		
  Erev	= -80	(mV)		
  Deadtime = 1	(ms)		
  GMAX     = 1  (mho)		
  DELAY = 0                     
}
INCLUDE "sns.inc"