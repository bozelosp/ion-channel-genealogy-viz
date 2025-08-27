NEURON {
	POINT_PROCESS mGluR_DynSyn
	USEION ca WRITE ica	
	RANGE tau_rise, tau_decay
	RANGE U1, tau_rec, tau_fac
	RANGE i, g, e, mg, inon, ica, ca_ratio
	NONSPECIFIC_CURRENT i
    }
    
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
    }    
    
    PARAMETER {
  	tau_rise 	= 20.0  (ms)  
	tau_decay 	= 1000.0 (ms) 
	U1        	= 1.0   (1)   
	tau_rec   	= 0.1   (ms)  
	tau_fac   	= 0.1   (ms)  
	e         	= 0.0   (mV)  
	mg			= 1.0   (mM)  
	ca_ratio  	= 0.1   (1)   
    }
    
    
ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	factor
	ica (nA)
	inon
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
	i = g*mgblock(v)*(v-e)
	ica=ca_ratio*i
	inon=(1-ca_ratio)*i
}

DERIVATIVE state{
	A' = -A/tau_rise
	B' = -B/tau_decay
}

FUNCTION mgblock(v(mV)) {
	
	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

NET_RECEIVE (weight, Pv, P, Use, t0 (ms)){
	INITIAL{
		P=1
		Use=0
		t0=t
		}	
	Use = Use * exp(-(t-t0)/tau_fac)
	Use = Use + U1*(1-Use) 
	P = 1-(1- P) * exp(-(t-t0)/tau_rec)
	Pv= Use * P
	P = P - Use * P
	
	t0=t
	
	A=A + weight*factor*Pv
	B=B + weight*factor*Pv
}