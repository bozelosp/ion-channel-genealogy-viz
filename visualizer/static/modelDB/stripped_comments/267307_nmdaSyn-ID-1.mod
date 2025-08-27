NEURON {
	POINT_PROCESS Exp2SynNMDA
	RANGE tau1, tau2, e, i, mgBlock, alpha_vspom, v0_block, eta, extMgConc, Kd, gamma, sh, mg_unblock_model
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>     
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	alpha_vspom = -0.062 (/mV) 
	                                           
	v0_block = 10 (mV)
	eta = 0.2801 (1)
	extMgConc = 1 (mM) 

	
	Kd = 9.888 (mM)
	gamma = 0.09137 (/mV)
	sh = 2.222 (mV)

	mg_unblock_model = 1 (1)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	mgBlock
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	mgBlock = vspom(v)
	i = g*mgBlock*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}

FUNCTION vspom (v(mV))( ){
	if (mg_unblock_model == 1) {
	   vspom = 1. / (1. + eta * extMgConc * exp(alpha_vspom * (v - v0_block))) 
	}
        else if (mg_unblock_model == 2) {
	   vspom = 1. / (1. + (extMgConc / 3.57) * exp(-0.062 * v))                
	}
	else if (mg_unblock_model == 3) {
  	   vspom = 1. / (1. + (extMgConc / Kd) * exp(gamma * (sh - v)))            
	}
}