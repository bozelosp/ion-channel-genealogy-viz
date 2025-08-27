TITLE Hyperpolarization activated inward current

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

NEURON {
    SUFFIX h
    NONSPECIFIC_CURRENT i
    RANGE gmax, g,q10
    THREADSAFE
}

PARAMETER {
    gmax  = 0.0002 (S/cm2)
	Eh = -25 (mV)
	q10 = 1 
}

ASSIGNED { 
    v (mV)
    g (S/cm2)
    i (mA/cm2)
	celsius (degC)
	m_inf
	tau_m
	qt
}

STATE {
    m
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	
    g = gmax*m
    i = g*(v-Eh)
}

INITIAL {
	settables(v)
	
    m = m_inf
    qt = q10^((celsius-6.3 (degC))/10 (degC))
}

DERIVATIVE states {
	settables(v)
	
    m' = (m_inf-m)/tau_m
}

UNITSOFF

PROCEDURE settables(v (mV)) {
	m_inf = 1/(1+exp((v+75)/12.5))
	tau_m = 3800/qt

}

UNITSON

