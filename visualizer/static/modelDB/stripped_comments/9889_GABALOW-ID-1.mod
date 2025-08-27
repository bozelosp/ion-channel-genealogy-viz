NEURON {  POINT_PROCESS GABALOW }









PARAMETER {
  Cmax	= 1	(mM)		
  Cdur	= 1.0	(ms)		
  Alpha	= 0.53	(/ms mM)	
  Beta	= 0.18	(/ms)		
  Erev	= -80	(mV)		
  Prethresh = 0 			
  Deadtime = 1	(ms)		
  gmax		(umho)		
}
INCLUDE "synq.inc"