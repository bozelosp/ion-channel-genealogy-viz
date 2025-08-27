INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iar
	USEION h READ eh WRITE ih VALENCE 1
	USEION ca READ cai
        RANGE ghbar, h_inf, tau_s, m, shift, o1, o2, p0, p1
	GLOBAL k2, cac, k4, Pc, nca, nexp, ginc, taum
}

UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(msM)	= (ms mM)
}


PARAMETER {
	eh	= -20	(mV)
	celsius = 36	(degC)
	ghbar	= 2e-5 (mho/cm2)
	cac	= 0.006 (mM)		
	k2	= 0.0001 (1/ms)		
	Pc	= 0.01			
	k4	= 0.001	(1/ms)		
	nca	= 4			
	nexp	= 1			
	ginc	= 2			
	taum	= 20	(ms)		
	shift	= 0	(mV)		
}


STATE {
	c1	
	o1	
	o2	
	p0	
	p1	
}


ASSIGNED {
	v	(mV)
	cai	(mM)
	ih	(mA/cm2)
        gh	(mho/cm2)
	h_inf
	tau_s	(ms)
	alpha	(1/ms)
	beta	(1/ms)
	k1ca	(1/ms)
	k3p	(1/ms)
	m
	tadj
}


BREAKPOINT {
	SOLVE ihkin METHOD sparse

	m = o1 + ginc * o2

	ih = ghbar * m * (v - eh)
}

KINETIC ihkin {









	evaluate_fct(v,cai)

	~ c1 <-> o1		(alpha,beta)

	~ p0 <-> p1		(k1ca,k2)

	~ o1 <-> o2		(k3p,k4)

	CONSERVE p0 + p1 = 1
	CONSERVE c1 + o1 + o2 = 1
}





INITIAL {




        tadj = 3.0 ^ ((celsius-36 (degC) )/10 (degC) )

	evaluate_fct(v,cai)

	c1 = 1
	o1 = 0
	o2 = 0
	p0 = 1
	p1 = 0
}


UNITSOFF
PROCEDURE evaluate_fct(v (mV), cai (mM)) {

	h_inf = 1 / ( 1 + exp((v+75-shift)/5.5) )

	tau_s = (taum + 1000 / ( exp((v+71.5-shift)/14.2) + exp(-(v+89-shift)/11.6) ) ) / tadj

	alpha = h_inf / tau_s
	beta  = ( 1 - h_inf ) / tau_s

	k1ca = k2 * (cai/cac)^nca

	k3p = k4 * (p1/Pc)^nexp

}






PROCEDURE activation(v (mV), cai (mM)) { LOCAL cc

	evaluate_fct(v,cai)

	cc = 1 / (1 + (cac/cai)^nca ) 		

	m = 1 / ( 1 + beta/alpha + (cc/Pc)^nexp )

	m = ( 1 + ginc * (cc/Pc)^nexp ) * m
}

UNITSON