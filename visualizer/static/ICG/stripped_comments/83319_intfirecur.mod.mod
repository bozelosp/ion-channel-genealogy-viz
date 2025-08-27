NEURON {
	ARTIFICIAL_CELL IntFireCur
	RANGE tau, m, refrac, minf, gid
	
}

PARAMETER {
	tau = 5 (ms)   <1e-9,1e9>
	refrac = 5 (ms) <0,1e9>
	minf = 0
	gid = -1
}

ASSIGNED {
	m
	t0(ms)
	on
	cnt
}

INITIAL {
	m = 0
	t0 = t
	on = 1 
	net_send(firetime(), 1)
}

FUNCTION M() {
	if (on) {
		M = minf + (m - minf)*exp(-(t - t0)/tau)
	}else{
		M = 0
	}
}

NET_RECEIVE (w) {
	m = M()
	t0 = t
	if (flag == 1) { 
		if (m >= .999999) { 
			net_event(t)
			m = 0
			on = 0 
			net_send(refrac + 1e-6, 2) 
		}else{
			net_send(firetime(), 1)
		}
	}else if (flag == 2) { 
		on = 1
		m = 0
		net_send(firetime(), 1)
	}else{ 
		if (on) {
			m = m + w
			if (m >= 1) {
				net_move(t)
			}else{
				net_move(t+firetime())
			}
		}
	}	
}

FUNCTION firetime()(ms) { 
	if (minf > 1) {
		firetime = tau*log((minf-m)/(minf - 1))
	}else{
		firetime = 1e20
	}

}