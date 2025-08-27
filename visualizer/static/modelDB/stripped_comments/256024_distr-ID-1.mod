UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
}


NEURON {
	SUFFIX ds
	RANGE vmax
}

ASSIGNED {
	vmax (mV)
}

INITIAL {
	vmax=v
}


BREAKPOINT {
	if (v>vmax) {vmax=v}
}