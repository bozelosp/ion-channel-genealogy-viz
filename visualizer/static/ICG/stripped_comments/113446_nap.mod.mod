UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

NEURON {
	SUFFIX nap
	USEION k READ ko
	USEION na READ nai, nao, ena WRITE ina
	GLOBAL ina_p_h, tau_act, conc_half, helling
	RANGE gnabar, ina
}

UNITS {
	
	FARADAY		= 96485.309 (coul)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	celsius		(degC)
	gnabar=1e-6	(mho/cm2)	
	helling=-.765	(mM)		
	conc_half=7 	(mM)		
	ina_p_h = 25000	(ms)		
	tau_act = 6	(ms)
}

ASSIGNED { 
	ina	(mA/cm2)
	ena	(mV)
	v	(mV)	
	nai	(mM)		
	nao	(mM)		
	ko	(mM)
}

STATE { ma mb ha hb }		

BREAKPOINT {
	SOLVE nastate METHOD sparse
	
	ina = gnabar*ma*ma*ha*kdep(ko)*(v-ena) 
	
	
}

INITIAL {
	
	ma=m_inf(v)
	mb=1-ma
	ha=h_inf(v)
	hb=1-ha
	ina = gnabar*ma*ma*ha*kdep(ko)*(v-ena) 
}

LOCAL a1,a2,b1,b2

KINETIC nastate {
	a1 = m_a(v)
	a2 = m_b(v)
	b1 = h_a(v)
	b2 = h_b(v)

	~ mb <-> ma (a1, a2)
	~ hb <-> ha (b1, b2)
	CONSERVE ma + mb = 1
	CONSERVE ha + hb = 1
}

FUNCTION kdep(ko (mM)) {
	TABLE DEPEND conc_half, helling FROM 0 TO 150 WITH 150
	kdep=1+ 2/(1+exp((ko-conc_half)/helling))
}

FUNCTION m_a(v(mV)) {
	
	TABLE FROM -150 TO 150 WITH 200
	
	
	
	
	
	m_a = m_inf(v)/tau_act
}

FUNCTION m_inf(v) {
	TABLE FROM -150 TO 150 WITH 200
	m_inf=1/(1+(exp(-(v+39.7)/7.0)))
}

FUNCTION m_b(v(mV)) {
	
	TABLE FROM -150 TO 150 WITH 200
	
	m_b = (1-m_inf(v))/tau_act
}

FUNCTION h_a(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	h_a = (1/ina_p_h)*(0.128*exp((7-v-70)/18))
}

FUNCTION h_b(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	h_b = (1/ina_p_h)*4/(1+exp((30-v-70)/5))
}


FUNCTION h_inf(v) {
	TABLE FROM -150 TO 150 WITH 200
	h_inf=h_a(v)/(h_a(v)+h_b(v))
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*1*FARADAY*v/(R*(celsius+273.11247574)) 
	eco = co*efun(z)
	eci = ci*efun(-z)
	
	
	ghk = (.001)*1*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}