NEURON {
	SUFFIX na_olm
	USEION na READ ena WRITE ina
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gl, el, ina
	RANGE Ra_m, vhalfa_m, qa_m, Rb_m, vhalfb_m, qb_m
	RANGE Ra_h, vhalfa_h, qa_h, Rb_h, vhalfb_h, qb_h
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gnabar = 0 (S/cm2)
	gl = 0 (S/cm2)

	
	Ra_m = 0.1(/ms)
	vhalfa_m = 32(mV)
	qa_m = 10(mV)

	Rb_m = 4(/ms)
	vhalfb_m = 57(mV)
	qb_m = 18(mV)

	Ra_h = 0.07(/ms)
	vhalfa_h = 58(mV)
	qa_h = 20(mV)

	Rb_h = 1(/ms)
	vhalfb_h = 28(mV)
	qb_h = 10(mV)
}

ASSIGNED {
	v (mV)
	el (mV)
	ena (mV)
	ina (mA/cm2)
	il (mA/cm2)
	minf
	hinf
	taum (ms)
	tauh (ms)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*pow(m, 3)*h*(v - ena)
	il = gl*(v - el)
}

DERIVATIVE states {
	
	rates(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	} else {
		Exp = exp(x)
	}
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {
		vtrap = 1(/mV)*x/(Exp(x/y) - 1)
	}
}

PROCEDURE rates(v(mV)) {
	
	
	LOCAL alpha, beta

	alpha = Ra_m*vtrap(-(v + vhalfa_m), qa_m)
	beta  = Rb_m*Exp(-(v + vhalfb_m)/qb_m)
	taum  = 1/(alpha + beta)
	minf  = alpha*taum

	alpha = Ra_h*Exp(-(v + vhalfa_h)/qa_h)
	beta  = Rb_h/(1 + Exp(-(v + vhalfb_h)/qb_h))
	tauh  = 1/(alpha + beta)
	hinf  = alpha*tauh
}