UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
     (mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	dt (ms)
	v (mV)
     ehd=-45  			(mV)
	ghdbar=0.00015 		(S/cm2)	
	gamma_ih	
	seed		
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT Iqq
	RANGE Iqq,ghdbar
}

STATE {
	qq
}

ASSIGNED {
	Iqq (mA/cm2)
}

INITIAL {
	qq=alpha(v)/(beta(v)+alpha(v))
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Iqq = ghdbar*qq*(v-ehd)
}

FUNCTION alpha(v(mV)) {
	alpha = 0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)			
}

FUNCTION beta(v(mV)) {
	beta = 0.001*193*exp(v/33.1)			
}

DERIVATIVE state {     
	qq' = (1-qq)*alpha(v) - qq*beta(v)
}