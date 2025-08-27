NEURON {
	THREADSAFE
	SUFFIX h_ca
	NONSPECIFIC_CURRENT i
	USEION ca READ cai
    RANGE gmax, i, tau, m
	
	GLOBAL e, taumin, vhalf, s1, s2
}

UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(msM)	= (ms mM)
}


PARAMETER {
	e	= -35	(mV)
	gmax= 2e-4 	(mho/cm2)
	cac	= 0.006 (mM)		
	k2	= 0.0001 (1/ms)		
	Pc	= 0.01				
	k4	= 0.001	(1/ms)		
	nca	= 4					
	nexp = 1				
	ginc = 2				


	vhalf = -78 (mV)
	vh2 = -85	(mV)
    s1 = -13 	(mV)
    s2 = 14 	(mV)
    taumax = 1020 (ms)		
    taumin = 20 (ms)		
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
	i	(mA/cm2)

	h_inf
	tau		(ms)
	alpha	(1/ms)
	beta	(1/ms)
	k1ca	(1/ms)
	k3p	(1/ms)
	m
	
}


BREAKPOINT {
	SOLVE ihkin METHOD sparse

	m = o1 + ginc * o2

	i = gmax * m * (v - e)
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




       

	evaluate_fct(v,cai)

	c1 = 1-h_inf
	o1 = h_inf
	o2 = 0
	p0 = 1
	p1 = 0
}


UNITSOFF
PROCEDURE evaluate_fct(v (mV), cai (mM)) {

	h_inf = 1 / ( 1 + exp((vhalf-v)/s1) )


	tau = 2*taumax/( exp((vh2-v)/s2) + exp((vhalf-v)/s1) ) + taumin

	alpha = h_inf / tau
	beta  = ( 1 - h_inf ) / tau

	k1ca = k2 * (cai/cac)^nca

	k3p = k4 * (p1/Pc)^nexp

}





PROCEDURE activation(v (mV), cai (mM)) { LOCAL cc

	evaluate_fct(v,cai)

	cc = 1 / (1 + (cac/cai)^nca ) 		

	
	
	
	m = ( 1 + ginc * (cc/Pc)^nexp ) / ( 1 + beta/alpha )
}

UNITSON