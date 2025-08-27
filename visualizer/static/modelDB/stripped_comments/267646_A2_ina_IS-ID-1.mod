NEURON {
	SUFFIX ina_AII_IS
	USEION na READ ena WRITE ina
	RANGE gnabar
	GLOBAL minf, hinf, mtau, htau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar = 0.2 (mho/cm2) 
	v_init = -62 (mV)
	mtau = 0.01 (ms)
	htau = 0.5 (ms)
	vhalfm_na = -48 (mV)
	vhalfh_na = -49.5 (mV)
	km_na = 5 (mV)
	kh_na = 2 (mV)
}

STATE {
	m h
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	ena (mV)
	minf hinf
	
}

INITIAL {
	m = 1/(1 + exp(-(v_init - vhalfm_na)/km_na))
	h = 1/(1 + exp((v_init - vhalfh_na)/kh_na))
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*m*m*m*h*(v - ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(v(mV)) {
	minf = 1/(1 + exp(-(v - vhalfm_na)/km_na))
	hinf = 1/(1 + exp((v - vhalfh_na)/kh_na))
}