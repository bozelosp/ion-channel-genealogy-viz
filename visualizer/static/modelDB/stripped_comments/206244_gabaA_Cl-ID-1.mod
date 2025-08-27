NEURON {
	POINT_PROCESS gaba

	USEION cl READ ecl WRITE icl VALENCE -1

	NONSPECIFIC_CURRENT ihco3
	RANGE onset, tau, gmax 
	RANGE P, HCO3e, HCO3i, i

	RANGE icl, ihco3, ehco3, e
}

UNITS {	
	(mA)    = (milliamp)
	(nA)    = (nanoamp)
	(mV)    = (millivolt)
	(uS)  = (micromho)
	(mM)    = (milli/liter)
	F = (faraday) (coulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {
	onset=0 	(ms)
	tau=.5 (ms)	<1e-3,1e6>
	gmax=0 	(uS)	<0,1e9>

	HCO3e   = 26	(mM)	
	HCO3i   = 16	(mM)	
	P       = 0.18		

	celsius = 34    (degC)
}

ASSIGNED {
	v	(mV)		

	icl	(nA)		
	ihco3	(nA)		
	i	(nA)		
				
	g 	(uS)		
				

	ecl	(mV)		
	ehco3	(mV)		
		
	e	(mV)		
}


INITIAL { 
	
	ehco3 = log(HCO3i/HCO3e)*(1000)*(celsius + 273.15)*R/F
	e = P*ehco3 + (1-P)*ecl
}

BREAKPOINT {
	
	if (gmax) { at_time(onset) }
	
	g = gmax * alphasyn(t - onset)
	
	icl = (1-P)*g*(v-ecl)

	ihco3 = P*g*(v-ehco3)
	i = icl + ihco3
	e = P*ehco3 + (1-P)*ecl

}

FUNCTION alphasyn(x) {
	if (x < 0 || x > (10 * tau)) {
		alphasyn = 0
	}else{
		alphasyn = x/tau * exp(-(x - tau)/tau)
	}
}