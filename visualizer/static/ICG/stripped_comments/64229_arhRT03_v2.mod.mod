NEURON {
    SUFFIX arhRT03 }
NEURON {
    
    NONSPECIFIC_CURRENT i
    GLOBAL eh
}

ASSIGNED {
	
	
	eh (mV)
}

PARAMETER {
  
  gmax 		= 0.4  (S/cm2)
  vrest           = 0    (mV)

  mvhalf 	= 75.
  mkconst 	= 5.5
  exptemp 	= 37.
  mq10		= 1
  mexp 		= 1

  hvhalf 	= 0
  hkconst 	= 0
  hq10		= 1
  hexp 		= 0
} 

INCLUDE "custom_code/inc_files/64229_boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { 
    settau = 1./(exp(-14.6 -0.086*v)+exp(-1.87 + 0.07*v))
  } else {
    settau = 1
  }
}

PROCEDURE iassign () {
	i=g*(v-eh)
}