INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX kca
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gmax, g, i, o_rate, c_rate, o
	GLOBAL cadep, maxc_rate, cainit
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(mA)	= (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gmax = 0.015	(mho/cm2)
	
	cadep = 2
	maxc_rate = 0.1 (1/ms)
	
	v		(mV)
	o_rate = 0	(1/ms)
	c_rate = 0	(1/ms)
	cainit = 5e-5	(mM)
}

ASSIGNED {
        ek      (mV)
        cai     (mM) 
	ik	(mA/cm2) 
	i	(mA/cm2)
        g       (mho/cm^2)
	o		
}

STATE {
	c		
}

BREAKPOINT {

	rates(cai)
	SOLVE state METHOD cnexp
	o = 1-c
        g = gmax*o
	i = g*(v-ek) ik=i 
}

DERIVATIVE state {

	c' = (c_rate*o) - (o_rate*c)
}

PROCEDURE rates(cai) {	

	TABLE o_rate, c_rate DEPEND cainit, cadep, maxc_rate FROM cainit TO 0.1 WITH 200

	o_rate = (cai - cainit) * cadep
	if (o_rate > 0) { c_rate = 1/o_rate
		if (c_rate > maxc_rate) { c_rate = maxc_rate }
			}
	else { c_rate = maxc_rate  }
}


INITIAL {
    rates(cai)
    o = o_rate/(o_rate+c_rate)
    c = 1-o
}