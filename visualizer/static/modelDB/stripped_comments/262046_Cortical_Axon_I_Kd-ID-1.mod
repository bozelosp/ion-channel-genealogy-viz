UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cortical_axon_i_kd
	USEION k WRITE ik				
	RANGE g_Kd, i_Kd				
}

PARAMETER {
	ek = -90 (mV)
	i_Kd = 0.0 (mA/cm2)				
	g_Kd = 0.6e-3 (S/cm2)
	tau_m = 1 (ms)
	tau_h = 1500 (ms)
	V_half_m = -43 (mV)
	V_half_h = -67 (mV)
	q_m = 8
	q_h = 7.3
	Q_s = 3.209						
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	m_inf
	h_inf
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_Kd*m*h*(v - ek)
	i_Kd = ik 						
}

UNITSOFF

INITIAL {
	settables(v)
	m = m_inf
	h = h_inf
}

DERIVATIVE states {
	settables(v)
	m' = (m_inf-m)/tau_m
	h' = (h_inf-h)/tau_h
}

PROCEDURE settables(v) {
	TABLE m_inf, h_inf FROM -100 TO 100 WITH 400
	
	m_inf = 1-(1/(1+exp((v-V_half_m)/q_m)))
	h_inf = 1-(1/(1+exp((v-V_half_h)/q_h)))
	
}

UNITSON