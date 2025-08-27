NEURON {
	POINT_PROCESS ExpSynSTDP
	RANGE tau, e, i, d, p, dtau, ptau, verbose, learning, LR, maxWeight, minWeight
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = 0	(mV)
	d = 0 <0,1>
	p = 0 
	dtau = 34 (ms) 
	ptau = 17 (ms) 
	verbose = 0
	learning = 1
	LR = 0.0001
	maxWeight = 1
	minWeight = 0
}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
}

STATE {
	g (uS)
}

INITIAL {
	g=0
	tpost = -1e9
	net_send(0, 1)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau
}

NET_RECEIVE(w (uS), tpre (ms)) {
	INITIAL { tpre = -1e9 }
	if (flag == 0) { 
		g = g + w
		if(learning) {
			if (w<=maxWeight){
				w = w+LR*p*0.15
				if (w>maxWeight){
        			w = maxWeight
        		}
			}
		}
		tpre = t
	}else if (flag == 2) { 
		tpost = t
		FOR_NETCONS(w1, tp) { 
        	if(learning) {
        		if (w1<=maxWeight){
	        		w1 = w1+LR*p*0.15
	        		if (w1>maxWeight){
	        			w1 = maxWeight
	        		}
	        		if(verbose) {
	        			printf("pot
	        		}
	        	}
        	}
		}
	} else { 
		WATCH (v > -20) 2
	}
}