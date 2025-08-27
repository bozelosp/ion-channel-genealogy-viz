INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kfast 
	THREADSAFE
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar, vshift, timefactor_n, ik
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
	gbar = 1   	(S/cm2)	
	v 		(mV)
	vshift = 0	(mV)
								
	tha  = 25	(mV)		
	qa   = 9	(mV)		
	
	Ra   = 0.02	(/ms)		
	Rb   = 0.002	(/ms)		

	dt		(ms)
	celsius		(degC)
	temp = 23	(degC)		
	q10  = 2.3			

	vmin = -120	(mV)
	vmax = 1000	(mV)
	
	timefactor_n = 1
} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(S/cm2)
	ek		(mV)
	ninf
	ntau (ms)	
	tadj
}
 

STATE { n }

INITIAL { 
	trates(v-vshift)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = tadj*gbar*n
	ik = gk * (v - ek)
} 



DERIVATIVE  states {   
        trates(v-vshift)      
        n' =  (ninf-n)/(timefactor_n*ntau)
}

PROCEDURE trates(v) {  
                      
        
        TABLE ninf, ntau
	DEPEND  celsius, temp, Ra, Rb, tha, qa
	
	FROM vmin TO vmax WITH 1600

	rates(v)





}


PROCEDURE rates(v) {  
                      

        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))

        tadj = q10^((celsius - temp)/10)
        ntau = 1/tadj/(a+b)
	if (ntau<1e-7) {
		ntau=1e-7
		
	}
	ninf = a/(a+b)
}