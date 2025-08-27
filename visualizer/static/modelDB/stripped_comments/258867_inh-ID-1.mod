NEURON {
	POINT_PROCESS inhSyn
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE g
	RANGE vgat,sst,npy,pv,xEff,V50,slope_factor
	RANGE isOn
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=-70	(mV)
	vgat=0
	sst=0
	npy=0
	pv=0
	xEff=-1
	isOn=0
    V50=-52 (mV)
    slope_factor=3
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
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
	g = rect(v)*(B - A)*isOn
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}


FUNCTION rect (v(mV))( ){
	rect= 1+(0.25-1)/ ( 1. + exp (( v - V50 )/slope_factor ) ) 
}