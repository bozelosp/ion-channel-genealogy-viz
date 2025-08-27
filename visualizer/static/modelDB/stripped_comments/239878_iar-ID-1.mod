INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
THREADSAFE
	SUFFIX iar
	USEION other WRITE iother VALENCE 1
        RANGE ghbar,  iother
	GLOBAL h_inf, tauh, erev, stp,  shift
}


UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(msM)	= (ms mM)
}


PARAMETER {
      v               (mV)
	erev	= -44	(mV)
	celsius = 36	(degC)
	ghbar	= 1.1e-5 (mho/cm2) 
	shift   =  0    (mV)

	
	
      stp     = 10	
	a0 = 96
	a1 = 250
	a2 = 30.7
	a3 = 78.8
	a4 = 5.78
}


STATE {
        h
}


ASSIGNED {
	i	(mA/cm2)
	iother 	(mA/cm2)
	h_inf
	tauh	(ms)
	tadj
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	iother = ghbar * h * (v - erev)
}

DERIVATIVE state  {
	evaluate_fct(v)
      h' = (h_inf - h) / tauh
}
UNITSOFF

INITIAL {
	evaluate_fct(v)
      h = h_inf
}


PROCEDURE evaluate_fct(v (mV)) {
	h_inf = 1 / ( 1 + exp((v+shift+a0)/stp) )
	tauh = exp((v+shift+a1)/a2) / ( 1 + exp((v+shift+a3)/a4))
}



UNITSON