INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX h
	NONSPECIFIC_CURRENT i
    RANGE gbar
    GLOBAL inf,tau
    GLOBAL eh
}

PARAMETER {
	gbar=8.0			(pS/um2) 	
									
	
	vhalf=-81			(mV)		
    a0=0.00057			(/ms)		
	zeta=7				(1)
    ab=0.4   			(1)
	qten=4.5			(1)			
	v					(mV)
    dt					(ms)
    temp = 33			(degC)		
	gas = 8.315			(J/degC)	
	farad = 9.648e4		(coul)		
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(S) = (siemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(J) = (joules)
}

ASSIGNED {
	eh (mV)
	i (mA/cm2)
    inf
    tau (ms)
	celsius (degC)		
}

STATE { 
		hh
}

INITIAL {
	rate(v)
	hh=inf
}

FUNCTION alpha(v(mV)) (1/ms) {
	alpha = a0*exp((0.001)*(-ab)*zeta*(v-vhalf)*farad/(gas*(273.16+celsius))) 
}

FUNCTION beta(v(mV)) (1/ms) {
	beta = a0*exp((0.001)*(1-ab)*zeta*(v-vhalf)*farad/(gas*(273.16+celsius)))
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    i = (0.0001)*gbar*hh*(v-eh)
}

DERIVATIVE state {
    rate(v)
    hh' = (inf-hh)/tau
}

PROCEDURE rate(v (mV)) { 
    LOCAL a, b, q10
    q10 = qten^((celsius-temp)/10(degC))
    a = q10*alpha(v)
    b = q10*beta(v)
    inf = a/(a+b)
    tau = 1/(a+b)
    if (tau<2) {tau=2}
}