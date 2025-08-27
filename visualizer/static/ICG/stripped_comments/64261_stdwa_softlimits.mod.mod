NEURON {
	POINT_PROCESS StdwaSoft
	RANGE interval, tlast_pre, tlast_post, M, P
	RANGE deltaw, wmax, aLTP, aLTD, wprune
	GLOBAL tauLTP, tauLTD, on
	POINTER wsyn
}

ASSIGNED {
	interval	(ms)	
	tlast_pre	(ms)	
	tlast_post	(ms)	
	M			
	P			
	deltaw			
	wsyn			
}

INITIAL {
	interval = 0
	tlast_pre = 0
	tlast_post = 0
	M = 0
	P = 0
	deltaw = 0
}

PARAMETER {
	tauLTP  = 20	(ms)    
	tauLTD  = 20	(ms)    
	wmax    = 1		
	aLTP    = 0.001		
	aLTD    = 0.00106	
	on	= 1		
	wprune  = 0             
}

NET_RECEIVE (w) {
	if (w >= 0) {				
		P = P*exp((tlast_pre-t)/tauLTP) + aLTP
		interval = tlast_post - t	
		tlast_pre = t
		deltaw = wsyn * M * exp(interval/tauLTD)
	} else {				
		M = M*exp((tlast_post-t)/tauLTD) - aLTD
		interval = t - tlast_pre	
		tlast_post = t
		deltaw = (wmax-wsyn) * P * exp(-interval/tauLTP)
	}
	if (on) {
		if (wsyn > wprune) {
		  wsyn = wsyn + deltaw
		} else {
		  wsyn = 0
		}
	}
}