INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE n, gkca, ikca, gbar
	RANGE ninf, ntau
	GLOBAL Ra, Rb, caix
	GLOBAL q10, temp, tadj, vmin, vmax
	THREADSAFE
}

UNITS {
	
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	
	gbar = 10   	(pS/um2)	
	v 				(mV)
	cai				(mM)
	caix = 4
									
	Ra   = 0.05		(/ms)		
	Rb   = 0.1		(/ms)		

	dt				(ms)
	celsius			(degC)
	temp = 23		(degC)		
	q10  = 2.3					

	vmin = -120		(mV)
	vmax = 100		(mV)
} 


ASSIGNED {
	
	a				(/ms)
	b				(/ms)
	ik 				(mA/cm2)
	ikca 			(mA/cm2)
	gkca			(pS/um2)
	ek				(mV)
	ninf
	ntau 			(ms)
	tadj
}
 

STATE { n }

INITIAL { 
	rates(cai)
	n = ninf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	gkca =gbar*n
	ikca = (1e-4) * gkca * (v - ek)
	ik = ikca
} 

LOCAL nexp

DERIVATIVE states {   
	
	rates(cai)
	n' = (ninf-n)/ntau
}

PROCEDURE rates(cai (mM)) {  
		
		a = Ra * cai^caix
        b = Rb

        tadj = q10^((celsius - temp)/10)

        ntau = 1/tadj/(a+b)
		ninf = a/(a+b)

 


}