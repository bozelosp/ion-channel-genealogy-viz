INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iar
	
        NONSPECIFIC_CURRENT i
        RANGE ghbar
	GLOBAL h_inf, tauh, stp, shift, eh
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
        eh (mV)
	i	(mA/cm2)
	
	h_inf
	tauh	(ms)
	tadj
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	i = ghbar * h * (v - eh)
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