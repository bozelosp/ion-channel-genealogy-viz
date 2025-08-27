INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX ito
	USEION k READ ek, ko, ki WRITE ik
	USEION ca READ cai
	RANGE g, ik, minf, mtau
}

PARAMETER {
	Kmca = 0.0005	(mM)
	Kmto = 10	(mM)
	g = 0.28	(uS/mM) 

}

STATE { 
	m 
}

ASSIGNED {
	v (mV)
	celsius (degC) 
	ik (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)
	ko (mM)
	ki (mM)
	cai (mM)
}

LOCAL RT
INITIAL {
	RT = (1000)*R*(273.15+celsius)
	rate(v)
	m = minf
}

BREAKPOINT { 
	SOLVE states METHOD derivimplicit



	ik = (1e-06)* m * g/S * (0.2 + (ko / (Kmto + ko))) * (cai / (Kmca + cai)) * (v + 10) / (1 - exp(-0.2(/mV) * (v + 10))) * (ki*exp(0.5 * v *F/RT) - ko*exp(-0.5 * v*F/RT))
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau
}

FUNCTION alp(v(mV)) (/ms) { 
	alp = (0.001)* 0.033(/s)*exp(-v / 17(mV))
}

FUNCTION bet(v(mV)) (/ms) { 
	bet = (0.001)* 33(/s) / (1 + exp(-(v + 10 (mV))/8(mV)))
}


PROCEDURE rate(v(mV)) { LOCAL a, b
TABLE minf, mtau FROM -100 TO 100 WITH 200
	a = alp(v)  b = bet(v) 
	mtau = 1/(a + b)
	minf = a * mtau
}