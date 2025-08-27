NEURON {
	SUFFIX cat
	USEION ca READ eca WRITE ica
	RANGE G, g, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	g = 50.97e-6 (siemens/cm2) <0,1e9>
}

ASSIGNED {
	v (mV)
	ek (mV)
	eca
	i
	ica (mA/cm2)
	G (siemens/cm2)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	G = g*m*h
	i = G*(v-eca)
	ica = i
}

INITIAL {
	m = minf(v)
	h = hinf(v)
}

DERIVATIVE states {
	m' = (minf(v)-m)/taum(v)
	h' = (hinf(v)-h)/tauh(v)
}


FUNCTION minf(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	
	minf = 1/(1+exp((v+20)/-1.898))
	
	
}

FUNCTION taum(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	
	
	taum = 18.51 - 3.388/(exp((v-6.53)/9.736)+exp((v+12.39)/-2.525))
}

FUNCTION hinf(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	
	hinf = 1/(1+exp((v+55.27)/6.11))    
}

FUNCTION tauh(v(mV)) {
	TABLE FROM -150 TO 150 WITH 500
	
	tauh = 94.16 - 45.27/(1 + exp((v+10)/-10))
	
	

}