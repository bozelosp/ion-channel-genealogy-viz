NEURON {
       SUFFIX bk
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gbar, gk,  ik, minf, taum, hinf, tauh, zinf, tauz
       RANGE zhalf,	ctauz
	   RANGE htau_factor,mtau_factor,ztau_factor,zinf_factor 
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
	
	cvm = 28.9 (mV)
	ckm = 6.2 (mV)

	ctm = 0.000505 (s)
	cvtm1 = 86.4 (mV)
	cktm1 = -10.1 (mV)
	cvtm2 = -33.3 (mV)
	cktm2 = 10 (mV)


	ch = 0.085
	cvh = 32 (mV)
	ckh = -5.8 (mV)
	cth = 0.0019 (s)
	cvth1 = 48.5 (mV)
	ckth1 = -5.2 (mV)
	cvth2 = -54.2 (mV)
	ckth2 = 12.9 (mV)

}

PARAMETER {
	v (mV)
	celsius (degC)
	ctauz = 1 (ms)
	mtau_factor = 1
	htau_factor = 1
	ztau_factor = 1

	gbar = 40 (pS/um2)

	ek (mV)
	cai (mM)

	zhalf = 0.01 (mM)
}

ASSIGNED {
	ik (mA/cm2)
    gk (pS/um2)   
	minf
	taum (ms)
	hinf
	tauh (ms)
	zinf
    tauz (ms)
}

STATE {
	m   FROM 0 TO 1
	z   FROM 0 TO 1
	h   FROM 0 TO 1
}

INITIAL {
	rates(v)
	m = minf
	z = zinf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    gk = gbar *  m^3 * z^2* h  
	ik = (1e-4)* gk * (v - ek)
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
	taum = (1e3) * ( ctm + 1 (s) / ( exp(-(v+cvtm1)/cktm1) + exp(-(v+cvtm2)/cktm2) ) ) / mtau_factor
	
	zinf =  1 /(1 + zhalf/cai)
    tauz = ctauz/ztau_factor

	hinf = ch + (1-ch) / ( 1+exp(-(v+cvh)/ckh) )
	tauh = (1e3) * ( cth + 1 (s) / ( exp(-(v+cvth1)/ckth1) + exp(-(v+cvth2)/ckth2) ) ) / htau_factor
}