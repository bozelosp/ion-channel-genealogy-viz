NEURON { SUFFIX nafRT03 }
NEURON {  USEION na READ ena WRITE ina }
ASSIGNED { ina }

PARAMETER {
	erev 		= 50    (mV)
	gmax 		= 0.4  (S/cm2)
        vrest           = 0    (mV)

	mvhalf 		= 34.5 
	mkconst 	= -10
	exptemp 	= 37
	mq10		= 1
	mexp 		= 3

	hvhalf 		= 62.9 
	hkconst 	= 10.7
	hq10		= 1
	hexp 		= 1
        ena
} 

INCLUDE "boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { 
    if (v<-26.5) { 
      settau = .025 + .14*exp((v+26.5)/10.) 
    } else {
      settau = .02 + .145*exp((-v-26.5)/10.) 
    }
  } else {
    settau = 0.15 + 1.15/(1.+exp((v+37.)/15.)) 
  }
}

PROCEDURE iassign () { i = g*(v-ena) ina=i }