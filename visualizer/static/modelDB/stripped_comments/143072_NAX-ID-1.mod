INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nax
	USEION na READ ena WRITE  ina		
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT inap
	
	NONSPECIFIC_CURRENT il
	
	
	

	RANGE gnapbar, gnabar, gl, ena, el
	RANGE mp_inf, m_inf, h_inf
	RANGE tau_mp, tau_m, tau_h
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gnapbar = 0.01	(mho/cm2)
	gnabar	= .2	(mho/cm2) 
	
	gl	= 0.007 (mho/cm2) 
	ena     = 55.0  (mV)
	
	
	el	= -71 (mV)		 
	
	dt              (ms)
	v               (mV)
	vtraub=-80
	ampA = 0.01
	ampB = 27		 
	ampC = 10.2	
	bmpA = 0.00025		
	bmpB = 34		
	bmpC = 10	
	amA = 1.76 	
	amB = 21.4	
	amC = 9.3	
	bmA = 0.13	
	bmB = 18.7 	
	bmC = 9.16
	ahA = 0.062	
	ahB = 114.0
	ahC = 11.0
	bhA = 1.7	
	bhB = 31.8
	bhC = 13.4
	
	
	
	
	
	
}

STATE {
	
	mp m h
}

ASSIGNED {
	
	inap    (mA/cm2)
	ina	(mA/cm2)
	
	il      (mA/cm2)
	mp_inf
	m_inf
	h_inf
	
	tau_mp
	tau_m
	tau_h
	
	q10_1
	q10_2
	q10_3
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	inap = gnapbar * mp*mp*mp * (v - ena)
	ina = gnabar * m*m*m*h * (v - ena)




	
	il   = gl * (v - el)
}

DERIVATIVE states {   
       evaluate_fct(v)
	mp'= (mp_inf - mp) / tau_mp
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	
}

UNITSOFF

INITIAL {




	
	
	

	evaluate_fct(v)
	mp = mp_inf
	m = m_inf
	h = h_inf
	
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2



























	a = vtrap1(v)
	b = vtrap2(v)
	tau_mp = 1 / (a + b)
	mp_inf = a / (a + b)

	a = vtrap6(v)
	b = vtrap7(v)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)

	a = vtrap8(v)
	b = bhA / (1 + Exp(-(v+bhB)/bhC))
	tau_h = 1 / (a + b)
	h_inf = a / (a + b)

	

	
	
	
	








}









FUNCTION vtrap1(x) {
	if (fabs((x+ampB)/ampC) < 1e-6) {
		vtrap1 = ampA*ampC
	}else{
		vtrap1 = (ampA*(x+ampB)) / (1 - Exp(-(x+ampB)/ampC))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+bmpB)/bmpC) < 1e-6) {
		vtrap2 = -bmpA*bmpC
	}else{
		vtrap2 = (bmpA*(-(x+bmpB))) / (1 - Exp((x+bmpB)/bmpC))
	}
}

FUNCTION vtrap6(x) {
	if (fabs((x+amB)/amC) < 1e-6) {
		vtrap6 = amA*amC
	}else{
		vtrap6 = (amA*(x+amB)) / (1 - Exp(-(x+amB)/amC))
	}
}

FUNCTION vtrap7(x) {
	if (fabs((x+bmB)/bmC) < 1e-6) {
		vtrap7 = -bmA*bmC
	}else{
		vtrap7 = (bmA*(-(x+bmB))) / (1 - Exp((x+bmB)/bmC))
	}
}

FUNCTION vtrap8(x) {
	if (fabs((x+ahB)/ahC) < 1e-6) {
		vtrap8 = -ahA*ahC
	}else{
		vtrap8 = (ahA*(-(x+ahB))) / (1 - Exp((x+ahB)/ahC)) 
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON