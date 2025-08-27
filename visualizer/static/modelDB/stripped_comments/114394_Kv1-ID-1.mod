NEURON {
	SUFFIX Kv1
	USEION k READ ek WRITE ik
	RANGE gk, gbar, ik
	GLOBAL ninf, taun
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
	q10 = 3

	ca = 0.10 (1/ms) 
	cva = 50 (mV)
	cka = -33.90877 (mV)

	cb = 0.12889 (1/ms)
     cvb = 50 (mV)
	ckb = 7.42101 (mV)         
}

PARAMETER {
	v (mV)
	celsius (degC)
	
	gbar = 0.011 (mho/cm2)   <0,1e9>
}


ASSIGNED {
 	ik (mA/cm2) 
	ek (mV)
	gk  (mho/cm2)
	ninf
	taun (ms)
	alphan (1/ms)
	betan (1/ms)
	qt
}

STATE { n }

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
      gk = gbar * n^4 
	ik = gk * (v - ek)
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