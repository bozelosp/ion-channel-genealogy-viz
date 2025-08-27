NEURON {
    SUFFIX rgi
    RANGE dc, sd, driver
    NONSPECIFIC_CURRENT i
}

UNITS {
    (mA) = (milliamp)
    (mA/cm2) = (nanoamp/cm2)
}


PARAMETER {
	dt	     (ms)
	dc	= 0. (mA/cm2)	
	sd	= 0. (mA/cm2)	
}

ASSIGNED {
    i (mA/cm2)              
    driver
}

INITIAL {
    i = dc
}

BREAKPOINT {
    i = dc - sd*driver
}