NEURON {
	SUFFIX Hor_Il
    NONSPECIFIC_CURRENT i
    RANGE i, e, gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (mS) = (millimho)
}

PARAMETER {
    gbar = 4.7170e-06 (mho/cm2) 
    e = -80 (mV)
}

ASSIGNED {
	i    (mA/cm2)
    v    (mV)
}

BREAKPOINT {
    i = gbar * (v - e)
}