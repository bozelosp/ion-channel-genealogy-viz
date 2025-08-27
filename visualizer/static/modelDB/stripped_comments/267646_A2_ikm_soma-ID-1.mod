NEURON {
	SUFFIX ikm_AII_soma
	USEION k READ ek WRITE ik
	RANGE gkmbar, ik
	GLOBAL minf, mtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkmbar = 0 (mho/cm2)
	v_init = -62 (mV)
	mtau = 50 (ms)
	vhalfm_km = -40 (mV)
	km_m = 4 (mV)
}

STATE {
	m
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ek (mV)
	minf
}

INITIAL {
	m = 1/(1 + exp(-(v_init - vhalfm_km)/km_m))
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkmbar*m*(v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
}

PROCEDURE rates(v (mV)) {
	minf = 1/(1 + exp(-(v - vhalfm_km)/km_m))
}