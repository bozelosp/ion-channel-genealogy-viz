INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Ikt3	USEION k READ ek WRITE ik
	RANGE n, gk, gbar, actshift, i
	GLOBAL ninf, ntau, ik
	GLOBAL Ra, Rb, tha, qa
	GLOBAL tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 41.8448   	(pS/um2)	
	v 		(mV)
								
	tha  = 40	(mV)		
	qa   = 16.7	(mV)		
	
	Ra   = 0.6	(/ms)		
	Rb   = 0.012	(/ms)		

	dt		(ms)
	celsius		(degC)
	

	vmin = -120	(mV)
	vmax = 60	(mV)
  actshift = 0 (mV)} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau (ms)	
	tadj
        i               (mA/cm2)
}
 

STATE { n }

INITIAL { 
	trates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n
	ik = (1e-4) * gk * (v +  98)
	i = ik

} 

LOCAL nexp

DERIVATIVE states {   
        trates(v)      
	n' = (ninf - n)/ntau
}

PROCEDURE trates(v) {  
                      
        TABLE ninf, ntau
	DEPEND celsius, Ra, Rb, tha, qa
	
	FROM vmin TO vmax WITH 199

	rates(v)

}


PROCEDURE rates(v) {  
                      

        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
        ntau = 1/(a+b)
	ninf = a*ntau
}