INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ittc
	USEION ca READ cai,cao WRITE ica
	GLOBAL q10m,q10h
	RANGE g, gmax, m_inf, tau_m, h_inf, tau_h, shift, i
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
	gmax	= 0.0022 (mho/cm2)
	q10m	= 3			
	q10h	= 3			
        exptemp = 24    (degC)
	shift	= 2 	(mV)		
	cai	= 2.4e-4 (mM)		
	cao	= 2	(mM)
}

STATE {
  m h
}

ASSIGNED {
	g	(mho/cm2)
	i	(mA/cm2)
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)			
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
        celsius
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	g = gmax * m * m * h
	i = g * (v-carev)
        ica = i
}

DERIVATIVE states {
	mh(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}


UNITSOFF
INITIAL {




	phi_m = q10m ^ ((celsius-exptemp)/10)
	phi_h = q10h ^ ((celsius-exptemp)/10)

	mh(v)
	h = h_inf
	m = m_inf
}

PROCEDURE mh (v(mV)) { LOCAL Vm

	Vm = v + shift

	m_inf = 1.0 / ( 1 + exp(-(Vm+57)/6.2) )
	h_inf = 1.0 / ( 1 + exp((Vm+81)/4.0) )


        tau_m = (1  /  (exp(-(Vm+129.6)/16.7) + exp((Vm+14.8)/18.2) ) + 0.612)/phi_m

        tau_h = (30.8+(211.4  +  exp((Vm+113.2)/5))/(1+exp((Vm+84)/3.2)))/phi_h
}

UNITSON