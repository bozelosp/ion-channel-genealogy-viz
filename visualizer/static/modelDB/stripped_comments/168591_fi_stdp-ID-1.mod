NEURON {
	POINT_PROCESS FastInhibSTDP

	
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE gmax
	RANGE mgid, ggid, srcgid

	
	RANGE interval, tlast_pre, tlast_post, M, P
	RANGE deltaw, wmax, aLTP, aLTD
	RANGE wsyn
	GLOBAL tauLTP, tauLTD, on
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	
	tau1=1 (ms) <1e-9,1e9>
	tau2 = 200 (ms) <1e-9,1e9>
	gmax = .003 (uS) 
	e = -80	(mV)

	
	tauLTP  = 20	(ms)    
	tauLTD  = 20	(ms)    
	wmax    = 1		
	aLTP    = 0.001		
	aLTD    = 0.00106	
	on	= 1		

	
	mgid = -1 
	ggid = -1 
	srcgid = -1 
}

ASSIGNED {
	
	v (mV)
	i (nA)
	g (uS)
	factor

	
	interval	(ms)	
	tlast_pre	(ms)	
	tlast_post	(ms)	
	M			
	P			
	deltaw			
	wsyn			

}

STATE {
	A
	B
}

INITIAL {
	
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	
	interval = 0
	tlast_pre = 0
	tlast_post = 0
	M = 0
	P = 0
	deltaw = 0
	wsyn = 0
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	g = (B - A)*gmax
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE (w) {
	if (w >= 0) {				
		P = P*exp((tlast_pre-t)/tauLTP) + aLTP
		interval = tlast_post - t	
		tlast_pre = t
		deltaw = wmax * M * exp(interval/tauLTD)
	} else {				
		M = M*exp((tlast_post-t)/tauLTD) - aLTD
		interval = t - tlast_pre	
		tlast_post = t
		deltaw = wmax * P * exp(-interval/tauLTP)
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
	A = A + wsyn*factor
	B = B + wsyn*factor
}