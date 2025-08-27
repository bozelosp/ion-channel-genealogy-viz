NEURON {
	SUFFIX ahp
	USEION ca READ cai
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar, q, tau_q, qinf, betaq_const
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 2167e-6	(S/cm2) < 0, 1e9 >
	Erev = -82 (mV)
	cai (mM) 
	betaq_const =  0.074	
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	qinf
	tau_q (ms)
}

STATE {	q }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * q*q
	i = g * (v - Erev)
}

INITIAL {
	
	q = alphaq(v)/(alphaq(v) + betaq(v))
}

DERIVATIVE states {
	rates(v)
	q' = (qinf - q)/tau_q
}

FUNCTION alphaq(Vm (mV)) (/ms) {
	UNITSOFF
	alphaq = 3.5e9 * (cai)^3
	UNITSON
}

FUNCTION betaq(Vm (mV)) (/ms) {
	UNITSOFF
	betaq = betaq_const
	UNITSON
}

FUNCTION tauq(Vm (mV)) (/ms) {
	UNITSOFF
	tauq = 1.0 / (alphaq(Vm) + betaq(Vm))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_q = tauq(Vm)
	qinf = alphaq(Vm) * tau_q      
}