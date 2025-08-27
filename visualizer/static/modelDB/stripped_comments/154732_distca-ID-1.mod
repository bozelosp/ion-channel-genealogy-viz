UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
}


NEURON {
	SUFFIX dca
	USEION ca READ cai	
        RANGE camax
}

PARAMETER {
	cai		(mM)
}

ASSIGNED {
	camax
}

INITIAL {
	camax=cai
}


BREAKPOINT {
	if (cai>camax) {camax=cai}
}