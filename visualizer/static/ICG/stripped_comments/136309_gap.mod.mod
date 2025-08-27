UNITS { 
	(mV) = (millivolt) 
	(nA) = (nanoamp) 
  (uS) = (microsiemens)
} 

NEURON { 
    POINT_PROCESS Gap
    POINTER vgap
	  RANGE g, i
	  NONSPECIFIC_CURRENT i
}

PARAMETER { 
	  g = 0.0 	   (uS)
} 

ASSIGNED {
    v (mV)
    vgap (mV)
    i (nA)
}

BREAKPOINT { 
	  i = g*(v - vgap) 
}