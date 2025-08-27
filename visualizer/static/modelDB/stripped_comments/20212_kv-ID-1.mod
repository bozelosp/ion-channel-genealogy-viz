INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kv
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar
	RANGE ninf, ntau
	GLOBAL Ra, Rb
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 5   	(pS/um2)	
	v 		(mV)
								
	tha  = 25	(mV)		
	qa   = 9	(mV)		
	
	Ra   = 0.02	(/ms)		
	Rb   = 0.002	(/ms)		

	dt		(ms)
	celsius		(degC)
	temp = 23	(degC)		
	q10  = 2.3			

	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau (ms)	
	tadj
}
 

STATE { n }

INITIAL { 
	trates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states
	gk = tadj*gbar*n
	ik = (1e-4) * gk * (v - ek)
} 

LOCAL nexp

PROCEDURE states() {   
        trates(v)      
        n = n + nexp*(ninf-n)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(v) {  
                      
        LOCAL tinc
        TABLE ninf, nexp
	DEPEND dt, celsius, temp, Ra, Rb, tha, qa
	
	FROM vmin TO vmax WITH 199

	rates(v)

        tadj = q10^((celsius - temp)/10)

        tinc = -dt * tadj
        nexp = 1 - exp(tinc/ntau)
}


PROCEDURE rates(v) {  
                      

        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
        ntau = 1/(a+b)
	ninf = a*ntau
}