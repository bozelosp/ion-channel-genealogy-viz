TITLE the delayed, rectifying potassium current

COMMENT
change it to Connor-Stevens model
ENDCOMMENT

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

NEURON {
    SUFFIX kd
    USEION k READ ko WRITE ik
    RANGE gmax,ik, q10
    THREADSAFE
}

CONSTANT {
	FARADAY = 96485.3399
	R = 8.31447215
}

PARAMETER {

    gmax = 0.003 (S/cm2)
	Ekd1 = -70 (mV)
	q10 = 1
	k_i = 140
}

ASSIGNED {
    v (mV)
    ko
    ek1 (mV)
    ek
    ik (mA/cm2)
	am
	bm
	m_inf
	tau_m
	celsius (degC)
	qt
}

STATE {
    m
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gmax*m*m*m*m*(v-Ekd1)
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
	am = (-.01*(v+50-4.3)/(exp(-(v+50-4.3)/10)-1))/2*3.8
	bm = .125*exp(-(v+60-4.3)/80)/2*3.8
	m_inf = am/(am+bm)
	tau_m = 1/(am+bm)*3.8/qt
}

UNITSON
