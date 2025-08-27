NEURON {
	SUFFIX Purk_lkg
	NONSPECIFIC_CURRENT i
	RANGE el, glbar, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	glbar = 0.002 (mho/cm2)
	el = -70 (mV)
}

ASSIGNED { i (mA/cm2) }

BREAKPOINT { i = glbar * (v - el ) }