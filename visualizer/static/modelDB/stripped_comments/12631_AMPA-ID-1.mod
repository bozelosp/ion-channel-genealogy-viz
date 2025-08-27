NEURON {  POINT_PROCESS AMPA }






PARAMETER {
  Cdur	= 1.1	(ms)		
  Alpha	= 10	(/ms mM)	
  Beta	= 0.5	(/ms)		
  Erev	= 0	(mV)		
  Deadtime = 2.5	(ms)		
  GMAX  = 1	(umho)		
  DELAY = 0                     
}
INCLUDE "sns.inc"