NEURON {
	POINT_PROCESS Exp2SynNmda
    NONSPECIFIC_CURRENT i
	RANGE tau1, tau2, e, i, mgBlock, extMgConc, alpha_vspom, v0_block 
	RANGE isOn

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	tau1= .1 (ms) <1e-9,1e9> 
	tau2 = 10 (ms) <1e-9,1e9> 
	e=0	(mV)
	alpha_vspom = -0.087 (/mV) 
	v0_block =  -3  (mV) 
	extMgConc = 1 
	isOn = 0
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
	i = isOn*g*mgBlock*(v - e)

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