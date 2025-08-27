NEURON {
	SUFFIX pcCaP
	USEION ca READ cai, cao WRITE ica
	RANGE pcabar, ica, ghk_value
	GLOBAL minf, taum
	GLOBAL monovalConc, monovalPerm
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	q10 = 2.2
	F = 9.6485e4 (coulombs)
	R = 8.3145 (joule/kelvin)

	cv = 19 (mV)
	ck = 5.5 (mV)
}

PARAMETER {
	v (mV)
	celsius (degC)

	cai (mM)
	cao (mM)

	pcabar = 6e-5 (cm/s)
	monovalConc = 140 (mM)
	monovalPerm = 0
}

ASSIGNED {
	qt
	ica (mA/cm2)
      minf 
	taum (ms)
	T (kelvin)
	E (volt)
	zeta
    ghk_value
}

STATE { m }

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	T = kelvinfkt( celsius )
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = (1e3) * pcabar * m * ghk(v, cai, cao, 2)
    ghk_value = ghk(v, cai, cao, 2)
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/taum
}

FUNCTION ghk( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) { 
	E = (1e-3) * v
      zeta = (z*F*E)/(R*T)	
	
	

	if ( fabs(1-exp(-zeta)) < 1e-6 ) {
	ghk = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)
	} else {
	ghk = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))
	}
}

PROCEDURE rates( v (mV) ) {
	minf = 1 / ( 1 + exp(-(v+cv)/ck) )
	taum = (1e3) * taumfkt(v)/qt
}

FUNCTION taumfkt( v (mV) ) (s) {
	UNITSOFF
	if ( v > -50 ) {
	taumfkt = 0.000191 + 0.00376 * exp(-((v+41.9)/27.8)^2)
	} else {
	taumfkt = 0.00026367 + 0.1278 * exp(0.10327*v)
	}
	UNITSON
}

FUNCTION kelvinfkt( t (degC) )  (kelvin) {
	UNITSOFF
	kelvinfkt = 273.19 + t
	UNITSON
}