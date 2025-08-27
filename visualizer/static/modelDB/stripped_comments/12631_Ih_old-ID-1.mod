INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iar
	USEION other WRITE iother VALENCE 1
	USEION ca READ cai
        RANGE ghbar, gh, i
	GLOBAL k2, cac, nexp, h_inf, tau_s, tau_f, controls, controlf
}

UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(msM)	= (ms mM)
}


PARAMETER {
	eh	= -43	(mV)
	celsius = 36	(degC)
	ghbar	= .0001	(mho/cm2)
	cac	= 1e-4	(mM)		
	k2	= 0.001	(1/ms)		
	nexp	= 2			
	controls = 1			
	controlf = 1			
}


STATE {
	s1
	s2
	f1
	f2
}


ASSIGNED {
	v	(mV)
	cai	(mM)
	i	(mA/cm2)
	iother	(mA/cm2)
        gh	(mho/cm2)
	h_inf
	tau_s	(ms)
	tau_f	(ms)
	alpha1	(1/ms)
	alpha2	(1/ms)
	beta1	(1/ms)
	beta2	(1/ms)
	kk	(1/ms)
	fderiv	(1/ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD runge

	if(controls == 0) {
		gh = ghbar * (f1+f2)
	} else if(controlf == 0) {
		gh = ghbar * (s1+s2)
	} else {
		gh = ghbar * (s1+s2) * (f1+f2)
	}
	
	i = gh * (v - eh)
	iother = i
}

DERIVATIVE states { LOCAL s0,f0
	evaluate_fct(v)

	s0 = 1 - s1 - s2
	f0 = 1 - f1 - f2

	kk = k2 * (5e-5/cac)^nexp

	fderiv = kk*s1 - k2*s2

	s1' = alpha1*s0 - beta1*s1 - fderiv
	s2' = fderiv

	fderiv = kk*f1 - k2*f2

	f1' = alpha2*f0 - beta2*f1 - fderiv
	f2' = fderiv
}

UNITSOFF
INITIAL {




        tadj = 3.0 ^ ((celsius-36)/10)
	evaluate_fct(v)
	kk = k2 * (cai/cac)^nexp
	s1 = alpha1*k2/(alpha1*kk + alpha1*k2 + beta1*k2)
	s2 = alpha1*kk/(alpha1*kk + alpha1*k2 + beta1*k2)
	f1 = alpha2*k2/(alpha2*kk + alpha2*k2 + beta2*k2)
	f2 = alpha2*kk/(alpha2*kk + alpha2*k2 + beta2*k2)
}


PROCEDURE evaluate_fct(v (mV)) {

	h_inf = 1 / ( 1 + exp((v+68.9)/6.5) )	
	tau_s = exp((v+183.6)/15.24) / tadj	
	tau_f = exp((v+158.6)/11.2) / ( 1 + exp((v+75)/5.5) ) / tadj

	alpha1 = controls * h_inf / tau_s
	beta1  = ( 1 - h_inf ) / tau_s
	alpha2 = controlf * h_inf / tau_f
	beta2  = ( 1 - h_inf ) / tau_f
}
UNITSON