NEURON{
	SUFFIX tonic
	USEION cl READ ecl WRITE icl VALENCE -1
	NONSPECIFIC_CURRENT ihco3
	RANGE  icl, leak, P, e,ihco3,ehco3 
	RANGE P, HCO3e, HCO3i, i
}

UNITS {
	(uS)  = (micromho)
	(nA)  = (nanoamp)
	(mV)  = (millivolt)
	(mM)    = (milli/liter)
	F 	  = (faraday) (coulombs)
	R     = (k-mole)  (joule/degC)
}


PARAMETER {
	leak    = 0.01 (siemens/cm2) <0, 1e9>
	P    = 0.18		
	celsius = 31    (degC)
	HCO3e   = 26	(mM)	
	HCO3i   = 16	(mM)	
}

ASSIGNED{
	icl (milliamp/cm2)
	i (milliamp/cm2)
	ihco3 (milliamp/cm2)
	v	(mV)		
	e	(mV)		
	ecl	(mV)		
    ehco3	(mV)		
}

INITIAL { 

	ehco3 = log(HCO3i/HCO3e)*(1000)*(celsius + 273.15)*R/F
	e = P*ehco3 + (1-P)*ecl
}


BREAKPOINT {
	icl = (1-P)*leak*(v-ecl)
	ihco3 = P*leak*(v-ehco3)
	i = icl + ihco3
	e = P*ehco3 + (1-P)*ecl
}