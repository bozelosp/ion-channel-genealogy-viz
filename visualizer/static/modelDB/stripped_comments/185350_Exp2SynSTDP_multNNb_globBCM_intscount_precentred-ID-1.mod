NEURON {
	POINT_PROCESS Exp2SynSTDP_multNNb_globBCM_intscount_precentred
	RANGE tau1, tau2, e, i, dtau, ptau
	RANGE wMax 	
	POINTER d, p	
	NONSPECIFIC_CURRENT i
	RANGE g, tpost, start
	
}

DEFINE EpspTimes 10000	

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	dtau = 36 (ms) 			
	ptau = 26 (ms) 			
	wMax = 0.01 (uS)		
	start = 30000 (ms)		
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	tpost (ms)
	presyntime[EpspTimes] (ms)	
	counter						
	total (uS)
	flagOLD
	d							
	p							
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	tpost = -100000

	FROM i=0 TO EpspTimes-1 {
		presyntime[i]=-1e9
	}
	counter = 0
	

	net_send(0, 1)
	flagOLD = 1
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}









NET_RECEIVE(w (uS), wE (uS), tpre (ms), X) {
	INITIAL { wE = w  tpre = -1e9	 X=0 }




	if (flag == 0) { 
		
		
		A = A + wE*factor
		B = B + wE*factor
		tpre = t
		counter=counter+1				
		presyntime[counter-1]=tpre		







		 
		if (t>start) {
			X = d*exp((tpost - t)/dtau)
			
			wE = wE*(1-X)
			if (wE>0) {} else {wE=0}
		}
		flagOLD = flag
	}else if (flag == 2) { 
		





		FOR_NETCONS(w1, wE1, tpres, X1) { 
		
		 
		if (flagOLD==flag) {} else {	
							if (t>start) {
								FROM i=0 TO counter-1 {
									X1 = p*exp((presyntime[i] - t)/ptau)
									wE1 = wE1*(1+X1)
								}
								
								if (wE1<wMax) {} else {wE1=wMax}	
								
								FROM i=0 TO counter-1 {
									presyntime[i]=-1e9
								}
								counter = 0	
							}
						    } 
		}
		tpost = t 
		
		flagOLD=flag
	} else { 

		








		WATCH (v > -37) 2 	
							
	}
}