NEURON {
	POINT_PROCESS tmgsyn
	RANGE e, i
	RANGE tau_1, tau_rec, tau_facil, U, u0
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	
	
	e = -90	(mV)
	
	
	tau_1 = 3 (ms) < 1e-9, 1e9 >
	
	
	tau_rec = 100 (ms) < 1e-9, 1e9 >
	
	
	tau_facil = 1000 (ms) < 0, 1e9 >
	
	
	
	
	U = 0.04 (1) < 0, 1 >
	
	u0 = 0 (1) < 0, 1 >
}

ASSIGNED {
	v (mV)
	i (nA)
	x
}

STATE {
	g (umho)
}

INITIAL {
	g=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau_1
}

NET_RECEIVE(weight (umho), y, z, u, tsyn (ms)) {
INITIAL {

	y = 0
	z = 0

	u = u0
	tsyn = t


}

	
	
	z = z*exp(-(t - tsyn)/tau_rec)
	z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )
	
	y = y*exp(-(t - tsyn)/tau_1)

	x = 1-y-z

	
	if (tau_facil > 0) {
		u = u*exp(-(t - tsyn)/tau_facil)
	} else {
		u = U
	}



	if (tau_facil > 0) {
		state_discontinuity(u, u + U*(1-u))
	}



	state_discontinuity(g, g + weight*x*u)
	state_discontinuity(y, y + x*u)

	tsyn = t


}