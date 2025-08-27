NEURON {
	SUFFIX KCa1
	USEION k READ ek WRITE ik
	USEION ca READ cai
	NONSPECIFIC_CURRENT i
	RANGE gbar, g,  ik, i, igate, nc
	GLOBAL minf, taum, hinf, tauh, zinf, tauz
	GLOBAL zhalf
	GLOBAL gateCurrent, gunit
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
	e0 = 1.60217646e-19 (coulombs)
	q10 = 2.7 (1)
	
	cvm = 28.9 (mV)
	ckm = 6.2 (mV)

	ctm = 0.000505 (s)
	cvtm1 = 86.4 (mV)
	cktm1 = -10.1 (mV)
	cvtm2 = -33.3 (mV)
	cktm2 = 10 (mV)

	ctauz = 1 (ms)

	ch = 0.085
	cvh = 32 (mV)
	ckh = -5.8 (mV)
	cth = 0.0019 (s)
	cvth1 = 48.5 (mV)
	ckth1 = -5.2 (mV)
	cvth2 = -54.2 (mV)
	ckth2 = 12.9 (mV)

	zm = 4.1023 (1)		
	zh = -4.3852 (1)		
}

PARAMETER {
	gateCurrent = 0 (1)	
	
	gbar = 0.007 (S/cm2)
	gunit = 182 (pS)		
	
	zhalf = 0.001 (mM)
}

ASSIGNED {
	celsius (degC)
	v	(mV)
	
	ik	(mA/cm2)
	i	(mA/cm2)
	igate (mA/cm2)
	g	(S/cm2)   
	
	ek	(mV)
	cai	(mM)

	nc	(1/cm2)		

	qt	(1)
	
	minf	(1)
	taum	(ms)
	hinf	(1)
	tauh	(ms)
	zinf	(1)
      tauz	(ms)
}

STATE {
	m   FROM 0 TO 1
	z   FROM 0 TO 1
	h   FROM 0 TO 1
}

INITIAL {
	nc = (1e12) * gbar / gunit
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	m = minf
	z = zinf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
      g = gbar * m^3 * z^2 * h      
	ik = g * (v - ek)
	igate = nc * (1e6) * e0 * ( 3 * zm * mgateFlip() + zh * hgateFlip() )

	if (gateCurrent != 0) { 
		i = igate
	}
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/taum
	z' = (zinf-z)/tauz
	h' = (hinf-h)/tauh
}

PROCEDURE rates( v (mV) ) {
	v = v + 5 (mV)
	minf = 1 / ( 1+exp(-(v+cvm)/ckm) )
	taum = (1e3) * ( ctm + 1 (s) / ( exp(-(v+cvtm1)/cktm1) + exp(-(v+cvtm2)/cktm2) ) ) / qt
	
	zinf = 1 /(1 + zhalf/cai)
      tauz = ctauz/qt

	hinf = ch + (1-ch) / ( 1+exp(-(v+cvh)/ckh) )
	tauh = (1e3) * ( cth + 1 (s) / ( exp(-(v+cvth1)/ckth1) + exp(-(v+cvth2)/ckth2) ) ) / qt
}

FUNCTION mgateFlip() (1/ms) {
	mgateFlip = (minf-m)/taum
}

FUNCTION hgateFlip() (1/ms) {
	hgateFlip = (hinf-h)/tauh
}