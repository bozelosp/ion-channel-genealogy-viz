TITLE Amacrine Cell Hodgkin Huxley Potassium Channel

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

NEURON {
	SUFFIX HHk
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL inf, tau
}

PARAMETER {
	gkbar = 0.0002					(mho/cm2)
	ek 								(mV)
}

STATE {
	n
}

ASSIGNED {
	v 								(mV)
	celsius 						(degC)
	ik 								(mA/cm2)
	inf
	tau 							(ms)
}

INITIAL {
	rate(v)
	n = inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar * n * n * n * n * (v - ek)
}

DERIVATIVE states {
	rate(v)
	n' = (inf - n)/tau
}	

UNITSOFF

FUNCTION alp(v(mV)) {
	LOCAL q10
	
	v = -v - 65
	q10 = 3
	alp = 0.01 * q10 * expM1(v + 10, 10)
}

FUNCTION bet(v(mV)) {
	LOCAL q10
	
	v = -v - 65
	q10 = 3^((celsius - 6.3)/10)
	bet = 0.125 * q10 * exp(v/80)
}

FUNCTION expM1(x, y) {
	if (fabs(x/y) < 1e-6) {
		expM1 = y * (1 - x/y/2)
	} 
	else {
		expM1 = x/(exp(x/y) - 1)
	}
}

PROCEDURE rate(v(mV)) {
	LOCAL a, b
	
	TABLE inf DEPEND celsius FROM -100 TO 100 WITH 200
	
	a = alp(v)/5
	b = bet(v)/5
	tau = 1/(a + b)
	inf = a/(a + b)
}

UNITSON