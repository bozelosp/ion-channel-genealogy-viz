NEURON {
	POINT_PROCESS Exp2SynNMDA
	USEION ca READ eca WRITE ica
	RANGE tau1, tau2, e, i, ica, mgBlock,theDrive,theEca
	NONSPECIFIC_CURRENT i,ioffset

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
	eca = 100 (mV)
	alpha_vspom = -0.062 (/mV) 
	                                           
	v0_block = 10 (mV)
	caComponent = 0.1 
	extMgConc = 1 (mM) 
}

ASSIGNED {
	v (mV)
	i (nA)
	ica (nA)
	ioffset (nA)
	g (uS)
	factor
	mgBlock
	theDrive (mV)
	theEca (mV)
	
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
	ica = caComponent*g*mgBlock*(v-eca)
	theDrive = v-eca 
	theEca = eca	
	ioffset = -ica
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
	vspom=1./(1.+0.2801*extMgConc*exp(alpha_vspom*(v-v0_block))) 
	
}