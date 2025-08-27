NEURON {
	POINT_PROCESS Exp2NMDA
	NONSPECIFIC_CURRENT i
	RANGE tau1, tau2, e, i, Mg, K0, delta, wf
	THREADSAFE
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
	(S)  = (siemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(J)  = (joules)
}

PARAMETER {

	tau1 = 8.8		(ms)	<1e-9,1e9>	
	tau2 = 500		(ms)

	celsius 		(degC)	

	Mg = 1			(mM)	
	K0 = 4.1		(mM)	
	delta = 0.8 	(1)		

	e = -0.7		(mV)	
}

CONSTANT {
	T = 273.16	(degC)
	F = 9.648e4	(coul)	
	R = 8.315	(J/degC)
	z = 2		(1)		
}

ASSIGNED {
	v		(mV)
	dt		(ms)
	i		(nA)
	factor
	wf
}

STATE {
	A
	B
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
	
	wf = 1
	Mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	i = (B - A)*Mgblock(v)*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight) {
	wf = weight*factor
	A = A + wf
	B = B + wf
}

FUNCTION Mgblock(v(mV)) {
	
	Mgblock = 1 / (1 + (Mg/K0)*exp((0.001)*(-z)*delta*F*v/R/(T+celsius)))
}