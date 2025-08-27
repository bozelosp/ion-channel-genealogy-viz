NEURON {
	SUFFIX Golgi_lkg
	NONSPECIFIC_CURRENT i
	RANGE el, glbar, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	glbar = 21e-6 (mho/cm2)
	celsius  (degC)
	el = -55 (mV)
}

ASSIGNED { i (mA/cm2) }

BREAKPOINT { i = glbar * (v - el ) }