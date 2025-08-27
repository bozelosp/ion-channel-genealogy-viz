NEURON {
	SUFFIX Kv3
	USEION k READ ek WRITE ik
	RANGE gbar, g, ik,vshift
	GLOBAL ninf, tau

}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(mS) = (millisiemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	e0 = 1.60217646e-19 (coulombs)
	q10 = 2.7

	ca = 0.22 (1/ms)
	cva = 16 (mV)
	cka = -26.5 (mV)
	cb = 0.22 (1/ms)
	cvb = 16 (mV)
	ckb = 26.5 (mV)
	
	zn = 1.9196 (1)		
}

PARAMETER {
	vshift = 0
	gbar = 0.005 (S/cm2)   <0,1e9>
}

ASSIGNED {
	celsius (degC)
	v (mV)
	
	ik (mA/cm2)
 
	ek (mV)
	g (S/cm2)
	qt (1)

	ninf (1)
	tau (ms)
	alpha (1/ms)
	beta (1/ms)
}

STATE { n }

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	rateConst(v)
	n = ninf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
      g = gbar * n^4 
	ik = g * (v - ek)

}

DERIVATIVE state {
	rateConst(v)
	n' = alpha * (1-n) - beta * n
}

PROCEDURE rateConst(v (mV)) {
	alpha = qt * alphaFkt(v)
	beta = qt * betaFkt(v)
	ninf = alpha / (alpha + beta) 
	tau = 1 / (alpha + beta)
}

FUNCTION alphaFkt(v (mV)) (1/ms) {
	alphaFkt = ca * exp(-(v+cva+vshift)/cka)
}

FUNCTION betaFkt(v (mV)) (1/ms) {
	betaFkt = cb * exp(-(v+cvb+vshift)/ckb)
}