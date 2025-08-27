TITLE Amacrine Cell Hodgkin Huxley Sodium Channel

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

NEURON {
	SUFFIX HHna
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
	GLOBAL infm, infh, taum, tauh
}

PARAMETER {
	gnabar = 0.0026 			(mho/cm2)
}

STATE {
	m
	h
}

ASSIGNED {
	v 							(mV)
	celsius 					(degC)
	ina 						(mA/cm2)
	ena 						(mV)
	infm
	infh
	taum
	tauh
}

INITIAL {
	rate(v)
	m = infm
	h = infh
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar * m * m * m * h * (v - ena)
}

DERIVATIVE states {
	rate(v)
	m' = (infm - m)/taum
	h' = (infh - h)/tauh
}

FUNCTION alp(v(mV), i) (/ms) {
	LOCAL q10
	
	v = -v - 65 (mV)
	q10 = 3
	
	if (i == 0) {
		alp = q10 * 0.1 (/ms) * expM1(v * 1 (/ms) + 25, 10)
	}
	else if (i == 1) {
		alp = q10 * 0.07 (/ms) * exp(v/20(mV))
	}
}

FUNCTION bet(v(mV), i) (/ms) {
	LOCAL q10
	
	v = -v - 65 (mV)
	q10 = 3^((celsius - 6.3 (degC))/10 (degC))
	
	if (i == 0) {
		bet = q10 * 4 (/ms) * exp(v/18 (mV))
	}
	else if (i == 1) {
		bet = q10 * 1 (/ms)/(exp(0.1(/mV) * v + 3) + 1)
	}
}

FUNCTION expM1 (x,y) {
	if (fabs(x/y) < 1e-6) {
		expM1 = y * (1 - x/y/2)
	}
	else {
		expM1 = x/(exp(x/y) - 1)
	}
}

PROCEDURE rate(v(mV)) {
	LOCAL a, b
	
	TABLE infm, infh, taum, tauh DEPEND celsius FROM -100 TO 100 WITH 200
	a = alp(v, 0)
	b = bet(v, 0)
	taum = 1/(a + b)
	infm = a/(a + b)
	
	a = alp(v, 1)/3.5
	b = bet(v, 1)/3.5
	tauh = 1/(a + b)
	infh = a/(a + b)
}