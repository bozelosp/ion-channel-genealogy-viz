INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kca_slow
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE n, gk, gbar, jikcas
	RANGE ninf, ntau
	GLOBAL Ra, Rb, caix
	GLOBAL q10, temp, tadj, vmin, vmax
	THREADSAFE Ra, Rb, caix, q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 4e-1   	(pS/um2)	
	v 		(mV)
	cai  		(mM)
	caix = 1
	cac = 2e-4
	cas = 1e-5
									
	Ra   = 3e-3	(/ms)		
	Rb   = 5e-5	(/ms)		

	dt		(ms)
	celsius		(degC)
	temp = 23	(degC)		
	q10  = 2.3			

	vmin = -120	(mV)
	vmax = 100	(mV)
	eK=  -95 (mV)
	jikcas 		(mA/cm2)
} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau 		(ms)	
	tadj
}
 

STATE { n }

INITIAL { 
	rates(cai)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = tadj*gbar*n
	ik = (1e-4) * gk * (v - eK)
	jikcas = ik
} 

LOCAL nexp

DERIVATIVE states {
        rates(cai)
        n' =  (ninf-n)/ntau
}

PROCEDURE rates(cai(mM)) {  
    a = Ra/(1+exp((cac-cai)/cas))
    b = Rb
    
    tadj = q10^((celsius - temp)/10)
    
    ntau = 1/tadj/(a+b)
    ninf = a/(a+b)
}