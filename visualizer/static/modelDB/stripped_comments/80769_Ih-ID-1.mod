NEURON {
	SUFFIX Ih
	NONSPECIFIC_CURRENT i
	RANGE i, ghbar, eh
	GLOBAL ninf, taun
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

CONSTANT {
	q10 = 3
	
	cvn = 90.1 (mV)
	ckn = -9.9 (mV)
		
	cct = 0.19 (s)
	cat = 0.72 (s)
	cvt = 81.5 (mV)
	ckt = 11.9 (mV)
}

PARAMETER {
	v (mV)
	celsius (degC)
	
	ghbar = 0.0002 (S/cm2)
	eh = -30 (mV)
}

ASSIGNED {
	i (mA/cm2)
	qt
	ninf
	taun (ms)
}

STATE { n }

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar * n * (v - eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/taun
}

PROCEDURE rates(v (mV)) {
	ninf = 1 / ( 1+exp(-(v+cvn)/ckn) )
	taun = (1e3) * ( cct + cat * exp(-((v+cvt)/ckt)^2) ) / qt
}