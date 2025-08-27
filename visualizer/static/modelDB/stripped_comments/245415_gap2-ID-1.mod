NEURON {
	
	POINT_PROCESS gap2
	NONSPECIFIC_CURRENT i
	RANGE i, g
	POINTER vgap	


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