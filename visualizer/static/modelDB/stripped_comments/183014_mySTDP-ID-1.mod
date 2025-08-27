NEURON {
	POINT_PROCESS mySTDP 
	RANGE tau, e, i, d, p, dtau, ptau 
	NONSPECIFIC_CURRENT i 
        GLOBAL verbose 
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = 0	(mV)
	d = 0.2 <0,1>
	p = 0.2 
	dtau = 34 (ms) 
	ptau = 17 (ms) 
        verbose = 0 
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

NET_RECEIVE(w (uS), A, tpre (ms)) { 
	INITIAL { A = 1  tpre = -1e9 }
	if (flag == 0) { 
          if(verbose) {printf("entry flag=%g t=%g w=%g A=%g tpre=%g tpost=%g\n", flag, t, w, A, tpre, tpost)}
		g = g + w*A
		tpre = t
		A = A * (1 - d*exp((tpost - t)/dtau))
	}else if (flag == 2) { 
          if(verbose) {printf("entry flag=%g t=%g tpost=%g\n", flag, t, tpost)}
		tpost = t 
		FOR_NETCONS(w1, A1, tp) { 
                  if(verbose) {printf("entry FOR_NETCONS w1=%g A1=%g tp=%g\n", w1, A1, tp)}
			A1 = A1*(1 + p*exp((tp - t)/ptau))
		}
	} else { 
          if(verbose) {printf("entry flag=%g t=%g\n", flag, t)}
		WATCH (v > -20) 2 
	}
}