NEURON {
	POINT_PROCESS GABAb_DynSyn	
	RANGE tau_rise, tau_decay
	RANGE U1, tau_rec, tau_fac
	RANGE i, g, e
	NONSPECIFIC_CURRENT i
}

PARAMETER {
	tau_rise  = 50   (ms)  
	tau_decay = 200   (ms)  
	U1        = 1.0   (1)   
	tau_rec   = 0.1   (ms)  
	tau_fac   = 0.1   (ms)  
	e         = -95   (mV)  
}
     

ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	factor
}

STATE {
	A	
	B	
}

INITIAL{
	LOCAL tp
	A = 0
	B = 0
	tp = (tau_rise*tau_decay)/(tau_decay-tau_rise)*log(tau_decay/tau_rise)
	factor = -exp(-tp/tau_rise)+exp(-tp/tau_decay)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B-A
	i = g*(v-e)
}

DERIVATIVE state{
	A' = -A/tau_rise
	B' = -B/tau_decay
}

NET_RECEIVE (weight, Pv, P, Use, t0 (ms)){
	INITIAL{
		P=1
		Use=0
		t0=t
	    }	
	Use = Use * exp(-(t-t0)/tau_fac)
	Use = Use + U1*(1-Use) 
	P   = 1-(1- P) * exp(-(t-t0)/tau_rec)
	Pv  = Use * P
	P   = P - Use * P
	
	t0 = t
	
	A = A + weight*factor*Pv
	B = B + weight*factor*Pv	
}