NEURON {
	POINT_PROCESS leak
	NONSPECIFIC_CURRENT i
	RANGE r, vrest
	}

UNITS {
	(mV) = (millivolt)
	(nA) = (nanoamp)
}

PARAMETER  {
	r = 15.7 (megaohm)
	vrest = -45 (mV)
	}

ASSIGNED {
	i	(nA)
	v	(mV)
	}


BREAKPOINT {
	i = (v - vrest)/r
	}