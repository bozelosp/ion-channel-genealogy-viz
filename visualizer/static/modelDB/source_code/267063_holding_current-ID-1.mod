

NEURON {
    POINT_PROCESS i_hold
	RANGE i_holding
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA)   = (nanoamp)

}

PARAMETER {
	i_holding=0
	}

ASSIGNED { 
	i (nA)  
	g_a	(umho) 
}


BREAKPOINT {
		i   = i_holding/1000*-1
}