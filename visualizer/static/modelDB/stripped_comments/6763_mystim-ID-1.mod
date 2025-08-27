INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS VClamp1
	RANGE vset, gain, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	vset = -65  (mV)
	gain = 1 (nA/mV)
	dt  (ms)
        v (mV)
}
ASSIGNED { i (nA) }

BREAKPOINT {
		i = gain*(vset-v)
}