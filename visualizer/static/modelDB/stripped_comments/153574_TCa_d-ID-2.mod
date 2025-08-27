INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX TCa_d
	USEION ca READ cai, cao WRITE ica
	RANGE gmax, shift 
	GLOBAL minf, mtau, hinf, htau 
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius	= 36	(degC) 

	gmax	= .0008	(mho/cm2)
	shift	= 0 	(mV)
	cai	
	cao	
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	minf
	mtau	(ms)
	hinf
	htau	(ms)
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD euler
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gmax * m*m*h * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}

UNITSOFF
INITIAL {






	cai	= 2.4e-4 (mM)		
	cao	= 2	(mM)    

	phi_m = 5.0 ^ ((celsius-24)/10)
	phi_h = 3.0 ^ ((celsius-24)/10)

	evaluate_fct(v)

	m = minf
	h = hinf
}

PROCEDURE evaluate_fct(v(mV)) { 




	minf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )
	hinf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )

	mtau = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m
	htau = ( 85 + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phi_h
}
UNITSON