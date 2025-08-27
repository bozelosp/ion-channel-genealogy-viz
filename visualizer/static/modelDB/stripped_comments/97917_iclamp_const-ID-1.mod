NEURON {
	POINT_PROCESS IClamp_const
	RANGE amp, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	amp (nA)
}
ASSIGNED { i (nA) }

INITIAL {
	i = amp
}

BREAKPOINT {
	
}