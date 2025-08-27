NEURON {
	ARTIFICIAL_CELL IntervalFire
	RANGE tau, m, invl, burst_start, burst_stop, burst_factor
	
	POINTER r
	RANGE noutput, ninput 
}

PARAMETER {
	tau = 5 (ms)   <1e-9,1e9>
	invl = 10 (ms) <1e-9,1e9> 
	burst_start = 0 (ms)
	burst_stop = 0 (ms)
	burst_factor = 1
}

ASSIGNED {
	m
	minf
	t0(ms)
	r
	tau1
	minf1
	ninput noutput
}

INITIAL {
	ninput = 0
	noutput = 0
	tau1 = 1/tau
	minf = 1/(1 - exp(-invl*tau1)) 
	minf1 = 1/(minf - 1)
	specify_invl() 
	m = 0
	t0 = t
	net_send(firetime(), 1)
}

FUNCTION M() {
	M = minf + (m - minf)*exp(-(t - t0)*tau1)
}

NET_RECEIVE (w) {
	m = M()
	t0 = t
	if (flag == 0) {
		ninput = ninput + 1
		m = m + w
		if (m > 1) {
			m = 0
			noutput = noutput + 1
			net_event(t)
		}
		net_move(t+firetime())
	}else{
		net_event(t)
		noutput = noutput + 1
		m = 0
		specify_invl()
		net_send(firetime(), 1)
	}
}

FUNCTION firetime()(ms) { 
	firetime = tau*log((minf-m)*minf1)

}

PROCEDURE specify_invl() {
VERBATIM {
	extern double nrn_random_pick(void*);
	if (!_p_r) {
		return 0.;
	}
	invl = nrn_random_pick((void*)_p_r);
	if (t >= burst_start && t <= burst_stop) {
		invl *= burst_factor;
	}
}
ENDVERBATIM
	minf = 1/(1 - exp(-invl*tau1)) 
	minf1 = 1/(minf - 1)
}

PROCEDURE set_rand() {
VERBATIM {
	extern void* nrn_random_arg(int);
	void** ppr;
	ppr = (void**)(&(_p_r));
	*ppr = nrn_random_arg(1);
}
ENDVERBATIM
}