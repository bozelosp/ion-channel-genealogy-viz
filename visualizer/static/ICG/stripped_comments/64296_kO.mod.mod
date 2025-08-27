NEURON {
	SUFFIX kO
	
	USEION k READ ek WRITE ik
        RANGE gkbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gkbar = 1.0 (mho/cm2)
	
}

ASSIGNED {
	ik (mA/cm2)
        ek (mV)
}

BREAKPOINT {
	ik = gkbar*(v - ek)
}