NEURON {
	SUFFIX nacx
	USEION na READ nai,nao WRITE ina
	USEION ca READ cao,cai WRITE ica
	RANGE ina,ica, icana
	RANGE knaca, dnaca, R, F, gbar
}

UNITS {
    (molar) = (1/liter)                     
    (mM)    = (millimolar)             	
	(uA) = (microamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(um) = (micron)

}

PARAMETER {
	knaca = 36e-9 
	dnaca = .0036 
	R = 8314 (mV)
	F = 96500  
	gbar = 316 (1/cm2) 

}

ASSIGNED {
	v (mV)
	nai (mM)
	nao (mM)
	cai (mM)
	cao (mM)
	ina (mA/cm2)
	ica (mA/cm2)
	celsius (degC)
	icana (mA/cm2)
	inaca (mA/cm2)
}


INITIAL {
	rates(v)
}

BREAKPOINT {
	rates(v)
	icana = inaca
	ina = 3*inaca		
	ica = -2*inaca

}

UNITSOFF
FUNCTION rates(Vm (mV)) {    
	LOCAL q10, T, dfcain, dfcaout, s
	T = 273 + celsius
	q10 = (2.2*(T-296.0)+(310.0-T))/14.0
	dfcain = nai^3*cao*exp(0.5*Vm*F/(R*T))
	dfcaout = nao^3*cai*exp(-0.5*Vm*F/(R*T))
	s = 1 + dnaca*(cai*nao^3 + cao*nai^3)
	inaca = gbar*q10*knaca*(dfcain-dfcaout)/s
}
UNITSON