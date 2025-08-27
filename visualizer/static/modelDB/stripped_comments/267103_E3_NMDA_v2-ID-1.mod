NEURON {
	POINT_PROCESS E3_NMDA_v2
	RANGE tau1, tau2, tau3, wtau2, wtau3, factor, e, i, scalar, v1
	RANGE g, B, C, E, open
	RANGE wf
	RANGE eta, gamma, mgblock
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	tau3 = 11 (ms) <1e-9,1e9>
	e = 0	(mV)
	factor = 10                 
	wtau2 = 0.5             
	

	Mg = 1	(mM)
	
	
	scalar = 1
}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
	g (uS)
	wf
	mgblock
	wtau3
	open
	v1
}

STATE {
	C (uS)
	B (uS)
	E (uS)

}

INITIAL {
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	if (tau2/tau3 > .9999) {
	    tau2 = .9999*tau3
    }
    wtau3 = 1-wtau2
	C = 0
	B = 0
	E = 0
	tpost = -1e9
	net_send(0, 1)
	mgblock = Mgblock(v)
	open = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = wtau2*B + wtau3*E - C
	mgblock = Mgblock(v)
	i = g*(v - e)*mgblock*scalar*0.17720458423318466 
	open = g * scalar
	v1 = v
}

DERIVATIVE state {
	C' = -C/tau1
	B' = -B/tau2
	E' = -E/tau3
      
}

NET_RECEIVE(w (uS)) {
	if (flag == 0) { 
		wf = w*factor
		C = C + wf
		B = B + wf
		E = E + wf
	}
}

FUNCTION Mgblock(v(mV)) {
	
	Mgblock = 1 / (1 + exp(-62 * v * 0.001(1/mV)) * Mg / 3.57)
}