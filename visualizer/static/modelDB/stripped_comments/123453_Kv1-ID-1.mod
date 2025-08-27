NEURON {
	SUFFIX Kv1
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT i
	RANGE g, gbar, ik, i , igate, nc
	GLOBAL ninf, taun
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

	ca = 0.12889 (1/ms)
	cva = 45 (mV)
	cka = -33.90877 (mV)

	cb = 0.12889 (1/ms)
      cvb = 45 (mV)
	ckb = 12.42101 (mV)

	zn = 2.7978 (1)		
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
	g  (S/cm2)
	nc (1/cm2)			
	
	ninf (1)
	taun (ms)
	alphan (1/ms)
	betan (1/ms)
	qt (1)
}

STATE { n }

INITIAL {
	nc = (1e12) * gbar / gunit
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
      g = gbar * n^4 
	ik = g * (v - ek)
	igate = nc * (1e6) * e0 * 4 * zn * ngateFlip()

	if (gateCurrent != 0) { 
		i = igate
	}
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/taun 
}

PROCEDURE rates(v (mV)) {
	alphan = alphanfkt(v)
	betan = betanfkt(v)
	ninf = alphan/(alphan+betan) 
	taun = 1/(qt*(alphan + betan))       
}

FUNCTION alphanfkt(v (mV)) (1/ms) {
	alphanfkt = ca * exp(-(v+cva)/cka) 
}

FUNCTION betanfkt(v (mV)) (1/ms) {
	betanfkt = cb * exp(-(v+cvb)/ckb)
}

FUNCTION ngateFlip() (1/ms) {
	ngateFlip = (ninf-n)/taun 
}