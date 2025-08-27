NEURON {
	SUFFIX dcnpas
	NONSPECIFIC_CURRENT i
	RANGE gbar, e, i
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gbar = .001	(S/cm2)	<0,1e9>
	e = -60	(mV) 
	    
	    
}

ASSIGNED {
    v (mV)  
    i (mA/cm2)
}

BREAKPOINT {
	i = gbar * (v - e)
}