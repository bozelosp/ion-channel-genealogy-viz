NEURON {

    POINT_PROCESS Gap
    POINTER vgap
    RANGE ggap, i
    NONSPECIFIC_CURRENT i
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER { ggap = 0 (nS) }
    
ASSIGNED {

    v    (mV)
    vgap (mV)
    i    (nA)
}
 
BREAKPOINT { 

	if (ggap>0) {i = (1e-3) * ggap * (v-vgap) }

}