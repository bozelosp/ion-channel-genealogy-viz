NEURON {
	SUFFIX L_Ca_inact
	
	USEION ca READ eca WRITE ica
	RANGE gcabar,ica,m_inf,m,h
	GLOBAL vca,theta_m,kappa_m,theta_h,kappa_h,eca
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gcabar  = 0.0003  (mho/cm2)
	eca		(mV)	

	dt		(ms)
	tau_m	= 20	(ms)
	v		(mV)
        vca=80		(mV)
	theta_m = -30   (mV)
	kappa_m = -6	(-mV)
	tau_h	= 1500	(ms)
	theta_h = 14   (mV)
	kappa_h = 4	(-mV)
}

STATE {
	m
	h
}

ASSIGNED {
	ica		(mA/cm2)
	m_inf
	h_inf
	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar * m *h* (v - eca)  
}

DERIVATIVE states {
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

	m_inf = 1 / (1 + (Exp((v - theta_m)/ kappa_m)))
	h_inf = 1 / (1 + (Exp((v - theta_h)/ kappa_h)))

}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(Exp(x/y)-1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}