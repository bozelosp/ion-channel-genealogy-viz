NEURON {
	SUFFIX ampa 
	POINT_PROCESS ampa
	RANGE tau1, tau2, e, i, g_max, g, k
	NONSPECIFIC_CURRENT i
	GLOBAL total, i2, g2
	EXTERNAL Area_canmda
}

UNITS {
	(nA) = (nanoamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1 = 0.5    (ms) <1e-9,1e9>   
	tau2 = 2.6   (ms) <1e-9,1e9>	
	g_max= 0.000010 (uS) <1e-9,1e9>	
	e    = 0       (mV) 
	
	Area  (cm2)
	k = 1e-06 (mA/nA)
}

ASSIGNED {
	v (mV)
	i (nA)
	factor
	total 
	g (uS)
	
	
	g2 (uS)	
	i2 (mA/cm2)
}

STATE {
	A 
	B 
}

INITIAL {
	LOCAL t_peak
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	t_peak = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-t_peak/tau1) + exp(-t_peak/tau2)
	factor = 1/factor

	Area = Area_canmda  
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = g_max*(B-A)
	i = g*(v-e)		
	
	g2=g
	i2=i*k/Area
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight) {
	state_discontinuity(A, weight*factor)
	state_discontinuity(B, weight*factor)
	total = total+weight
}