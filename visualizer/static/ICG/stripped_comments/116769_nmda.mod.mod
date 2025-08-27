NEURON {
	POINT_PROCESS nmda
	RANGE tau1, tau2, tau3, e, i, g_max, g, A, B, C	,k
	NONSPECIFIC_CURRENT i
	GLOBAL total,i2,g2 
	EXTERNAL Area_canmda
}

UNITS {
	(nA) = (nanoamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1 = 3.18     (ms) <1e-9,1e9>     
	tau2 = 57.14      (ms) <1e-9,1e9>	
	tau3 = 2000     (ms) <1e-9,1e9>	    
	
	g_max= 0.000045 (uS)			
	e    = 0 (mV)
	mg   = 1 (mM)

	Area	(cm2)
	k = 1e-06 (mA/nA)
}

ASSIGNED {
	v (mV)
	i (nA)
	factor
	total (uS)
	g (uS)
	
	g2 (uS)	
	i2 (mA/cm2)	
}

STATE {
	A (uS)
	B (uS)
	C (uS)
}

INITIAL {
	LOCAL t_peak
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	C = 0
	
	factor=0.8279	
	factor = 1/factor
	
	Area = Area_canmda 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	
	g = g_max*(B*0.8+C*0.2-A)
	i = g*(v - e)*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	
	
	g2=g			
	i2=i*k/Area		
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
	C' = -C/tau3
}

NET_RECEIVE(weight (uS)) {
	state_discontinuity(A, weight*factor)
	state_discontinuity(B, weight*factor)
	state_discontinuity(C, weight*factor)
	total = total+weight
	
}