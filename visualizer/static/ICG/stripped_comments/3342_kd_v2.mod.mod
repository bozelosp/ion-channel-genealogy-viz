INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kd
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar
	GLOBAL ninf, ntau, ik
	GLOBAL nna, nnc, nnq1, nnq2, ncen, nslp
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 74.338   	(pS/um2)	
	v 		(mV)

	ncen = -46.104	(mV)
	nslp = 1.4503	(mV)
	nna   = 5.8073	(/ms)		
	nnc   = -68.210	(mV)		
	nnq1  = 10.625	(mV)		
	nnq2  = 71.411	(mV)		

	dt		(ms)
	celsius		(degC)
	temp = 16	(degC)		
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
        SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n
	ik = (1e-4) * gk * (v - ek)
} 

LOCAL nexp

DERIVATIVE states {   
        trates(v)      
	n' = (ninf - n)/ntau
}

PROCEDURE trates(v) {  
                      
        TABLE ninf, ntau
	DEPEND celsius, temp, ncen, nslp, nna, nnc, nnq1, nnq2
	
	FROM vmin TO vmax WITH 199

	rates(v)

}


PROCEDURE rates(v) {  
                      

        ntau = nna / ( exp((v-nnc)/nnq1) + exp(-(v-nnc)/nnq2) )
	ninf = 1 / ( 1 + exp(-(v-ncen)/nslp) )
}