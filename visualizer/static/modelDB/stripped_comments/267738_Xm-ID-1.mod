NEURON {
	POINT_PROCESS Xm
	RANGE amp, cli
	USEION cl WRITE cli VALENCE 1
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	amp = 5	
}

ASSIGNED {
    cli (nA)
}

BREAKPOINT {
	cli = amp
}