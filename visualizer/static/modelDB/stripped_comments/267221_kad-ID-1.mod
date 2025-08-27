NEURON {
	SUFFIX kad
	USEION k READ ek WRITE ik
	RANGE gkabar, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
	gkabar = 0.018 (S/cm2)
}

ASSIGNED {
	v (mV)
	ek (mV)

	gka (S/cm2)
	ik (mA/cm2)
	ninf
	linf
	taun (ms)
	taul (ms)
}

STATE {
	n
	l
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gkabar*pow(n, 4)*l
	ik = gka*(v-ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/taun
	l' = (linf - l)/taul
}

INITIAL {
	rates(v)
	n = ninf
	l = linf
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION alpn(v (mV)) (/ms) {
	alpn = 0.01(/ms)*vtrap(-(v+34.4(mV)), 21(mV))
}

FUNCTION betn(v (mV)) (/ms) {
	betn = 0.01(/ms)*vtrap(v+34.4(mV), 21(mV))
}

FUNCTION alpl(v (mV)) (/ms) {
	alpl = -0.01(/ms)*vtrap(v+58(mV), 8.2(mV))
}

FUNCTION betl(v (mV)) (/ms) {
	betl = -0.01(/ms)*vtrap(-(v+58(mV)), 8.2(mV))
}

PROCEDURE rates(v (mV)) {
	
	ninf = alpn(v)/(alpn(v) + betn(v))
	taun = 0.2

	linf = alpl(v)/(alpl(v) + betl(v))
	if (v > -20) {
		taul = 5(ms) + 2.6(ms)*(v+20(mV))/10(mV)
	} else {
		taul = 5
	}
}