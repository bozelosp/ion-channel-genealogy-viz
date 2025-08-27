TITLE the transient, inward sodium current

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
    SUFFIX na
    USEION na READ nai,ena WRITE ina
    RANGE gmax,ina, q10,m,h
    THREADSAFE
}

PARAMETER {
    gmax  = 0.014 (S/cm2)
    q10 = 1
}

ASSIGNED {
    v (mV)
    ina (mA/cm2)
	
	nai
	ena
    Ena1 (mV)
	
	m_inf
	tau_m
	h_inf
	tau_h
	qt
    celsius (degC)
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gmax*m*m*m*h*(v-ena)
:    ina = gmax*m*m*m*h*(v-64)
}

INITIAL {
	settables(v)
    m = m_inf
    h = h_inf
    qt = q10^((celsius-6.3 (degC))/10 (degC))
}

DERIVATIVE states {
	settables(v)
    m' = (m_inf-m)/tau_m
    h' = (h_inf-h)/tau_h
}
UNITSOFF

PROCEDURE settables(v (mV)) {
	m_inf = 1/(1+exp(-(v+48-10)/8.5))
	tau_m = (0.132/cosh((v+27)/7.5)+0.003/(1+exp(-(v+27)/5)))/qt
	
	h_inf = 1/(1+exp((v+47)/6))
	tau_h = 10/cosh((v+42)/15)/qt

}

UNITSON

