NEURON {
	SUFFIX caL
	USEION ca READ cai, cao WRITE ica
	RANGE Pbar, P, i, mu
	GLOBAL minf, mtau
	
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {
	vh = -35	(mV)	
	ve = 6.1	(mV)	
	mtauconst = 0.1	(ms)	



	Pbar = 42	(nanometer/s)	<0,1e9>
	cao = 2		(mM)
	cai = 10e-6	(mM)
	mu = 1
}

ASSIGNED {
	celsius	(degC)
	v	(mV)
	i	(uA/cm2)	
	ica	(mA/cm2)
	minf	(1)
	mtau	(ms)
	P	(nanometer/s)
	zFRT	(1/volt)
	zVFRT	(1)
	ghk	(coulomb/liter)
	
}


STATE {
	m
}


INITIAL {
	zFRT = (1000)*2*FARADAY/(R*(273+celsius))
	rates(v)
	m = minf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	P = mu*Pbar*m
	zVFRT = (0.001)*zFRT*v
	ghk = 2*FARADAY*(cai - cao*exp(-zVFRT))*gtrap(zVFRT)
	i = (1e-4)*P*ghk
	ica = (1e-3)*i
}


DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
}




PROCEDURE rates(v(mV)) {
UNITSOFF
	
	mtau = mtauconst
	minf = 1/(1 + exp(-(v - vh)/ve))
}
UNITSON


FUNCTION gtrap(x) {
	if (fabs(x) < 1e-6) {
		gtrap = 1 + x/2
	} else {
		gtrap = x/(1 - exp(-x))
	}
}