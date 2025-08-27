NEURON {
	POINT_PROCESS Golgi_extraGABA
	NONSPECIFIC_CURRENT i
	RANGE g, i, e

}

PARAMETER {
	g = 216e-06 (microsiemens)
	e = -80 (millivolt)
}
ASSIGNED {
	i (nanoamp)
	v (millivolt)
}
BREAKPOINT {
	i = g*(v - e)
}