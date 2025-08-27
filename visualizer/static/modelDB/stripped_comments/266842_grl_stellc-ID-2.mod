NEURON {
	POINT_PROCESS StellC
	RANGE v_SC, g_ampa, tau_ampa, e_ampa, i_ampa
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	Cm = 4 (pF)
	epas = -56 (mV)
	Grest = 0.2 (nS)
	
	e_ampa = 0	(mV)
	tau_ampa = 0.64 (ms)
}

ASSIGNED {
	i_ampa (nA)
}

STATE {
	v_SC (mV)
	g_ampa (uS)
}

INITIAL {
	v_SC = epas
	g_ampa=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	if (v_SC>10 && v_SC<20) {
		v_SC = epas
		}
	i_ampa = g_ampa*(e_ampa - v_SC)
}

DERIVATIVE state {
	v_SC' = (i_ampa + Grest*(epas - v_SC))/Cm
	g_ampa' = -g_ampa/tau_ampa
}

NET_RECEIVE(weight (uS)) {
	if (weight>=1) {
		g_ampa = g_ampa + weight
	}
	
	if (weight>=0.01 && weight<=0.03 && v_SC>-40) {
		v_SC = 20
		net_event(t) 
		g_ampa=0
	}
}