NEURON {
	POINT_PROCESS Xm
	RANGE amp
	USEION cl WRITE cli VALENCE 1
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	amp = -8	
}

ASSIGNED {
    cli (nA)
}

BREAKPOINT {
	cli = amp
}