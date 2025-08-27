: A-type (fast) Potassium Channel for AII amacrine cell
: Parameters and equations adapted from H. Riecke 2014
: Intrinsic bursting of AII amacrine cells underlies oscillations in the rd1 mouse retina

TITLE A-type Potassium Channel

NEURON {
	SUFFIX ika_AII_soma
	USEION k READ ek WRITE ik
	RANGE gkabar, ik
	GLOBAL minf, mtau, hinf, h1tau, h2tau, c
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkabar = 0.004 (mho/cm2)
	v_init = -62 (mV)
	mtau = 1 (ms)
	vhalfm_ka = -10 (mV)
	vhalfh_ka = -40.5 (mV)
	vhalfc_ka = -45 (mV)
	km_a = 7 (mV)
	kh_a = 2 (mV)
	kc_a = 15 (mV)
	f = 0.83
}

STATE {
	m h1 h2
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ek (mV)
	c
	minf hinf
	h1tau h2tau
}

INITIAL {
	rates(v)
	m = 1/(1 + exp(-(v_init - vhalfm_ka)/km_a))
	h1 = f*(1/(1 + exp((v_init - vhalfh_ka)/kh_a))) + (1 - f)
	h2 = f*(1/(1 + exp((v_init - vhalfh_ka)/kh_a))) + (1 - f)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	c = 1/(1 + exp(-(v + 45)/15))
	ik = gkabar*m*(c*h1 + (1 - c)*h2)*(v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h1' = (hinf - h1)/h1tau
	h2' = (hinf - h2)/h2tau
}	

PROCEDURE rates(v (mV)) {
	h1tau = 25 - 1/(20*(1 + exp(-(v + 35)/6)))
	
	if (((v + 17)^2/4 + 26) < 100) {
		h2tau = (v + 17)^2/4 + 26
	}
	else {
		h2tau = 100
	}
	
	minf = 1/(1 + exp(-(v - vhalfm_ka)/km_a))
	hinf = f*(1/(1 + exp((v - vhalfh_ka)/kh_a))) + (1 - f)
}

