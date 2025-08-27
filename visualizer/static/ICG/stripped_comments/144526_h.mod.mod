UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
     (mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	dt (ms)
	v (mV)
     
	gbar = 0.00015 (mho/cm2)	
	
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT i
	RANGE gbar
        GLOBAL eh
}

STATE {
	h
}

ASSIGNED {
         eh (mV)
	i (mA/cm2)
}

INITIAL {
	h=alpha(v)/(beta(v)+alpha(v))
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = gbar*h*(v-eh)
}

FUNCTION alpha(v(mV)) {
	alpha = 0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)			
}

FUNCTION beta(v(mV)) {
	beta = 0.001*193*exp(v/33.1)			
}

DERIVATIVE state {     
	h' = (1-h)*alpha(v) - h*beta(v)
}