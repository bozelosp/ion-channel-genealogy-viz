TITLE Ih-current

COMMENT

Author: Stefan Hallermann 

Provides deterministic Ih-currents as described in Kole et al. (2006).
	
ENDCOMMENT



UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	dt (ms)
	v (mV)
      eh  			(mV) 		:ih-reversal potential			       
	ghdbar=0.00015 		(S/cm2)	:default Ih conductance; exponential distribution is set in Ri18init.hoc 
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT i
	RANGE ghdbar,gamma_ih,eh
}

STATE {
	qq
}

ASSIGNED {
	i (mA/cm2)
}

INITIAL {
	qq=alpha(v)/(beta(v)+alpha(v))
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = ghdbar*qq*(v-eh)
}

FUNCTION alpha(v(mV)) {
	alpha = 0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)			:parameters are estimated by direct fitting of HH model to activation time constants and voltage actication curve recorded at 34C
}

FUNCTION beta(v(mV)) {
	beta = 0.001*193*exp(v/33.1)			
}

DERIVATIVE state {
	qq' = (1-qq)*alpha(v) - qq*beta(v)
}
