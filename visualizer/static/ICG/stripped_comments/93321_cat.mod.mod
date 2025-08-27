INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cat
	USEION ca READ cai,cao WRITE ica
	
	RANGE m_inf, tau_m, h_inf, tau_h, shift, i, carev, gbar
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
        (S) = (siemens)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)


	gbar = 1.0 (S/cm2) 

	cai  (mM) 
                                        
	cao	= 3	(mM)            
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	i	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
        celsius (degC)
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gbar * m*m*m*h * (v-carev)
        i = ica
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 

	m_inf= 1.0 / (1+exp( -(v+27.1)/7.2 ))
	h_inf= 1.0 / (1+exp( (v+32.1)/5.5 ))

	tau_m =  21.7 - 21.3 / (1+exp( -(v+68.1)/20.5 ))
	tau_h =  105 - 89.8 / (1+exp( -(v+55)/16.9 ))

}
UNITSON