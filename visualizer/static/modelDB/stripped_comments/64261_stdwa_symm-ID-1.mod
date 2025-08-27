NEURON {
	POINT_PROCESS StdwaSymm
	RANGE interval, tlast_pre, tlast_post
	RANGE deltaw, wmax, f
	GLOBAL tau_a, tau_b, a, on
	POINTER wsyn
}

ASSIGNED {
	interval	(ms)	
	tlast_pre	(ms)	
	tlast_post	(ms)	
	f                       
	deltaw			
	wsyn			
	tas             (ms2)   
}

INITIAL {
	interval = 0
	tlast_pre = 0
	tlast_post = 0
	f = 0
	deltaw = 0
}

PARAMETER {
	tau_a   = 20 (ms)       
	tau_b   = 15 (ms) 	
	wmax    = 1		
	a       = 0.001		
	on	= 1		
}

NET_RECEIVE (w) {
	tas = tau_a * tau_a 

	if (w >= 0) {				
		interval = tlast_post - t
		tlast_pre = t
		f = (1 - interval*interval/tas) * exp(interval/tau_b)
		deltaw = wmax * a * f
	} else {				
		interval = t - tlast_pre
		tlast_post = t
		f = (1 - interval*interval/tas) * exp(-interval/tau_b)
		deltaw = wmax * a* f
	}
	if (on) {
		wsyn = wsyn + deltaw
		if (wsyn > wmax) {
			wsyn = wmax
		}
		if (wsyn < 0) {
			wsyn = 0
		}
	}
}