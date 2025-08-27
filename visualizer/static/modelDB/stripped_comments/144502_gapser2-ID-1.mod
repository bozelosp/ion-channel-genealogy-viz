NEURON {
	POINT_PROCESS gap2
	NONSPECIFIC_CURRENT i
	RANGE i, g,t1,t2,g1,g2
	POINTER vgap	
}
PARAMETER {
	v (millivolt)
	g = 0 (nanosiemens) 
	g1=0 (nanosiemens)
	g2=0 (nanosiemens)
	t1=0 (milliseconds)
	t2=0 (milliseconds)
}
ASSIGNED {
	i (nanoamp)
	vgap (millivolt)
}
BREAKPOINT {
	if(t < t1) {
	     g = g1
	     
	}
	if( (t > t1) && (t < t2)) {
	    g = g2
	    
	}
	if(t > t2) {
	    g = g1
	    
	}
	i = (v - vgap)*g*(0.001)
}