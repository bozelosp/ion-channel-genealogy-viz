NEURON {
	POINT_PROCESS NMDA_DynSyn
	USEION ca WRITE ica	
	USEION mg READ mgo VALENCE 2
	RANGE tau_rise, tau_decay
	RANGE U1, tau_rec, tau_fac, stp
	RANGE i, g, e, mgb, inon, ica, ca_ratio
	NONSPECIFIC_CURRENT inon
    }
    
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
    }    
    
    PARAMETER {
  	tau_rise  = 5.0   (ms)  
	tau_decay = 70.0  (ms)  
	U1        = 1.0   (1)   
	tau_rec   = 0.1   (ms)  
	tau_fac   = 0.1   (ms)  
	e         = 0.0   (mV)  
	mgo		  = 1.0   (mM)  
	ca_ratio  = 0.1   (1)   
	stp       = 1.0   (1)   
    }
    
    
ASSIGNED {
	v		(mV)
	i		(nA)
	g		(umho)
	factor	(1)
	ica		(nA)
	inon	(nA)
	mgb
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
	mgb = mgblock(v)
	i = g*mgb*(v-e)
	ica = ca_ratio*i
	inon = (1-ca_ratio)*i
	
}

DERIVATIVE state{
	A' = -A/tau_rise
	B' = -B/tau_decay
}

FUNCTION mgblock(v(mV)) {
	
	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mgo / 3.57 (mM)))
}

NET_RECEIVE (weight, Pv, P, Use, t0 (ms)){
	INITIAL{
		P=1
		Use=0
		t0=t
	}	

    if(stp){
        Use = Use * exp(-(t-t0)/tau_fac)
        Use = Use + U1*(1-Use)
        P   = 1-(1- P) * exp(-(t-t0)/tau_rec)
        Pv  = Use * P
        P   = P - Use * P

        t0 = t

        A = A + weight*factor*Pv
        B = B + weight*factor*Pv
    } else {
        A = A + weight*factor
        B = B + weight*factor
    }
}