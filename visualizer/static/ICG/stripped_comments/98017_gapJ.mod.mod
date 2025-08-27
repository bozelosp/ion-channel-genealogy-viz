NEURON {
	POINT_PROCESS gapJ
	NONSPECIFIC_CURRENT i
	RANGE g, i
	POINTER vgap
}

PARAMETER {
	g	(umho)
}

ASSIGNED {
	i 	(milliamp)
	v 	(millivolt)
	vgap 	(millivolt)
}

BREAKPOINT {
	i = g*(v - vgap)
}