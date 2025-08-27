NEURON {
	SUFFIX Kv4
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT i
	RANGE gbar, g, ik, i, igate, nc
	GLOBAL ninf, taun, hinf, tauh
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
	q10 = 2.7

	can = 0.15743 (1/ms)
	cvan = 57 (mV)
	ckan = -32.19976 (mV)
	cbn = 0.15743 (1/ms)
	cvbn = 57 (mV)
	ckbn = 37.51346 (mV)

	cah = 0.01342 (1/ms)
	cvah = 60 (mV)
	ckah = -7.86476 (mV)
	cbh = 0.04477 (1/ms)
	cvbh = 54 (mV)
	ckbh = 11.3615 (mV)

	zn = 1.4736 (1)		
	zh = -5.4726 (1)		
}

PARAMETER {
	gateCurrent = 0 (1)	
	
	gbar = 0.004 (S/cm2)   <0,1e9>
	gunit = 16 (pS)		
}

ASSIGNED {
	celsius (degC)
	v (mV)

	ik (mA/cm2)
	i (mA/cm2)
	igate (mA/cm2)
 
	ek (mV)
	g (S/cm2)
	nc (1/cm2)			
	qt (1)

	ninf (1)
	taun (ms)
	alphan (1/ms)
	betan (1/ms)

	hinf (1)
	tauh (ms)
	alphah (1/ms)
	betah (1/ms)        
}

STATE { n h }

INITIAL {
	nc = (1e12) * gbar / gunit
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
      g = gbar * n^4 * h 
	ik = g * (v - ek)
	igate = nc * (1e6) * e0 * ( 4 * zn * ngateFlip() + zh * hgateFlip() )

	if (gateCurrent != 0) { 
		i = igate
	}
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/taun
	h' = (hinf-h)/tauh 
}

PROCEDURE rates(v (mV)) {
	alphan = alphanfkt(v)
	betan = betanfkt(v)
	ninf = alphan / (alphan + betan) 
	taun = 1 / (qt*(alphan + betan))
	alphah = alphahfkt(v)
	betah = betahfkt(v)
	hinf = alphah / (alphah + betah)
	tauh = 1 / (qt*(alphah + betah))       
}

FUNCTION alphanfkt(v (mV)) (1/ms) {
	alphanfkt = can * exp(-(v+cvan)/ckan) 
}

FUNCTION betanfkt(v (mV)) (1/ms) {
	betanfkt = cbn * exp(-(v+cvbn)/ckbn)
}

FUNCTION alphahfkt(v (mV))  (1/ms) {
	alphahfkt = cah / (1+exp(-(v+cvah)/ckah))
}

FUNCTION betahfkt(v (mV))  (1/ms)  {
	betahfkt = cbh / (1+exp(-(v+cvbh)/ckbh))
}

FUNCTION ngateFlip() (1/ms) {
	ngateFlip = (ninf-n)/taun 
}

FUNCTION hgateFlip() (1/ms) {
	hgateFlip = (hinf-h)/tauh 
}