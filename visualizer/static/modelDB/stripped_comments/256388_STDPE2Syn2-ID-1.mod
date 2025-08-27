NEURON {
	POINT_PROCESS STDPE2
	RANGE tau1, tau2, e, i, d, p, dtau, ptau, thresh, wmax, wmin
	RANGE g, gbdel, gblen, gbint, gscale, factor,dM,dV,B,C
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
	e = 0	(mV)
     pi=3.14159
	wmax = 0.0015 (uS)
      
	wmin = 0.0005 (uS)	
	d = 8 
	p = 1.2 
	
       dM = -22   (ms)
       dV= 5    (ms)
	 ptau = 10 (ms) 

	
      thresh = -55 (mV)	
	gbdel = 50 (ms) <1e-9,1e9> 
	gbint = 125 (ms) <1e-9,1e9> 
	gblen = 125 (ms) <1e-9,1e9> 
	gscale = 0.4	
      
}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
	on
	g (uS)
	gs
	factor
}

STATE {
	C (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	C = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor    
	gs=1
	on=0	
	tpost = -1e9
	net_send(0, 1)
	net_send(gbdel, 3)	
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - C
     	i = g*gs*(v - e)
    
}

DERIVATIVE state {
	C' = -C/tau1
	B' = -B/tau2
}

NET_RECEIVE(w (uS), A, tpre (ms) ) {
	INITIAL { A = 0  tpre = -1e9 }
	if (flag == 0) { 
	

		C = C + (w + A)*factor
		B = B + (w + A)*factor
          	
		tpre = t
		if (on == 1) {
			A = A * (1-(d*exp(-((tpost-t)-dM)^2/(2*dV*dV))) /(sqrt(2*pi)*dV))
		}
	}else if (flag == 2 && on == 1) { 

		tpost = t
            FOR_NETCONS(w1, A1, tp) { 
			
			A1 = A1 + (wmax-w1-A1)*p*exp((tp - t)/ptau)
			}
	} else if (flag == 1) { 

		WATCH (v > thresh) 2
	}
	else if (flag == 3) { 
		if (on == 0) { 
			on = 1
			gs = gscale
			net_send(gblen, 3)
		}
		else { 
			on = 0
			gs = 1
			net_send(gbint, 3)
		}
	}
}