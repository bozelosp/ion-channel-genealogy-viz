NEURON {  POINT_PROCESS AMPA }






PARAMETER {
  Cmax	= 1	(mM)		
  Cdur	= 1.1	(ms)		
  Alpha	= 10	(/ms mM)	
  Beta	= 0.5	(/ms)		
  Erev	= 0	(mV)		
  Prethresh = 0 			
  Deadtime = 2.5	(ms)		
  gmax		(umho)		
}
INCLUDE "synq.inc"