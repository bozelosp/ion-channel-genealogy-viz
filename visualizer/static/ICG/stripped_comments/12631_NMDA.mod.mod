NEURON{ POINT_PROCESS NMDA
  RANGE B 
}

PARAMETER {
  mg    = 1.    (mM)            
  Cdur	= 1.	(ms)		
  Alpha	= 4.	(/ms mM)	
  Beta	= 0.0067 (/ms)		
  Erev	= 0.	(mV)		
  Deadtime = 1	(ms)		
  GMAX	= 1     (S/cm2)         
  DELAY = 0                     
}

ASSIGNED { B }
INCLUDE "sns.inc"

BREAKPOINT {
  rates(v)
  g = g * B
  i = i * B
}
PROCEDURE rates(v(mV)) {
  TABLE B
  DEPEND mg
  FROM -100 TO 80 WITH 180
  B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}