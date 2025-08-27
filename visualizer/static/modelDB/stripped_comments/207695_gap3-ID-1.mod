NEURON {
	THREADSAFE
	POINT_PROCESS gap3
	NONSPECIFIC_CURRENT i
	RANGE i, g
	RANGE vgap	


}
PARAMETER {
	v (millivolt)
	g = 0 (nanosiemens) 
}
ASSIGNED {
	i (nanoamp)
	vgap (millivolt)
}
BREAKPOINT {

	i = (v - vgap)*g*(0.001)
}