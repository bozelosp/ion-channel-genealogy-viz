INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX TC_iT_Des98
	USEION cat READ cati, cato WRITE icat
	RANGE pcabar, m_inf, tau_m, h_inf, tau_h, shift, actshift
	GLOBAL qm, qh, km, kh
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
	
	celsius 	(degC)  
	pcabar	=.2e-3	(cm/s)	
	shift	= 2 	(mV)	
	actshift = 0 	(mV)	
	cati	= 2.4e-4 (mM)	
	cato	= 2	(mM)
	qm      = 2.5		
	qh      = 2.5           
        km      = 6.2
        kh      = 4.0
}

STATE {
	m h
}

ASSIGNED {
	icat	(mA/cm2)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	icat = pcabar * m*m*h * ghk(v, cati, cato)
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}


UNITSOFF
INITIAL {
	phi_m = qm ^ ((celsius-24)/10)
	phi_h = qh ^ ((celsius-24)/10)

	evaluate_fct(v)

	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) {
















	m_inf = 1.0 / ( 1 + exp(-(v+shift+actshift+57)/km) )
	h_inf = 1.0 / ( 1 + exp((v+shift+81)/kh) )

	tau_m = ( 0.612 + 1.0 / ( exp(-(v+shift+actshift+132)/16.7) + exp((v+shift+actshift+16.8)/18.2) ) ) / phi_m
	if( (v+shift) < -80) {
		tau_h = exp((v+shift+467)/66.6) / phi_h
	} else {
		tau_h = ( 28 + exp(-(v+shift+22)/10.5) ) / phi_h
	}

	
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	
	
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
FUNCTION nongat(v,cati,cato) {	
	nongat = pcabar * ghk(v, cati, cato)
}
UNITSON