NEURON {
	SUFFIX fn
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnabar, gkbar, gkmodbar
	RANGE fastNashift, vha, vhb
	RANGE am, an, lamb
	RANGE minf, ninf, ntau
	RANGE totna, totk
}

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}


PARAMETER { 
	fastNashift = 0	(mV)	
	gnabar = 0.0	(mho/cm2)
	gkbar = 1.0	(mho/cm2)
	am = 0.055	(/mV)
	vhm = -33	(mV)
	an = 0.055	(/mV)
	vhn = -40	(mV)
	lamb = 0.1		(1)
	tea = 0
	washfac = -0.005 (/ms)
}

ASSIGNED { 
	v	(mV)
	ena	(mV)  
	ek	(mV)  
	ina 	(mA/cm2) 
	ik 	(mA/cm2) 
	totna 	(mA/cm2) 
	totk 		(mA/cm2) 
	minf	(1)
	ninf	(1)
	ntau	(ms)

}

STATE { 
        n 
	gkmodbar (mho/cm2)
}


INITIAL {
	rates(v)
	if (tea > 0) {
		gkmodbar = gkbar
	}
	n = ninf

}

BREAKPOINT {
	SOLVE states METHOD cnexp
	totna = gnabar * minf^3 * (1-n) * ( v - ena ) 
	totk = gkbar * n^4 * ( v - ek )
	ina = gnabar * minf^3 * (1-n) * ( v - ena ) 
	ik = gkbar * n^4 * ( v - ek )
	if (tea > 0) {
		
		ik = gkmodbar * n^4 * ( v - ek )
	}
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/ntau
	
	gkmodbar' = washfac * gkmodbar
}


UNITSOFF 

PROCEDURE rates(V (mV)) {
	minf = 1 / ( 1 + exp( -2 * am * ( V - fastNashift - vhm )) )
	ninf = 1 / ( 1 + exp( -2 * an * ( V - vhn )) )
	ntau = 1 / lamb / (exp( an*(V-vhn)) + exp (-an*(V-vhn)) )
}

UNITSON