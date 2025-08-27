INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX it2
	
	USEION ca READ cai,cao WRITE ica
        RANGE gcabar, g, shift1
	GLOBAL m_inf, tau_m, h_inf, tau_h, shift2, sm, sh, phi_m, phi_h, hx, mx,rat
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

	gcabar	= .024	(mho/cm2)
	shift1	= -1 	(mV)
        shift2  = -6    (mV) 
        sm      = 7.4
        sh      = 5.0
        hx      = 1.5
        mx      = 3.0
	cai	(mM)		
	cao	(mM)
	rat	= 1
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	g       (mho/cm2)
	
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	g = gcabar * m*m*h
	ica = g * ghk(v, cai, cao)
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
	VERBATIM
	cai = _ion_cai;
	cao = _ion_cao;
	ENDVERBATIM






	phi_m = mx ^ ((celsius-24)/10)
	phi_h = hx ^ ((celsius-24)/10)

	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 




	m_inf = 1.0 / ( 1 + exp(-(v+shift1+50)/sm) )
	h_inf = 1.0 / ( 1 + exp((v+shift2+78)/sh) )

	tau_m = (2+1.0/(exp((v+shift1+35)/10)+exp(-(v+shift1+100)/15)))/ phi_m
	tau_h = (24.22+1.0/(exp((v+55.56)/3.24)+exp(-(v+383.56)/51.26)))/phi_h
}

FUNCTION ghk(v(mV), Ci(mM), Co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = Co*efun(z)*rat
	eci = Ci*efun(-z)
	
	
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
UNITSON