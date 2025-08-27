INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kc
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE n, gk, gbar
	RANGE ninf, ntau, hinf, htau
	GLOBAL Ra, Rb, Rc, Rd, caix
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {

	gbar = 1.0   	(pS/um2)	
	v 		(mV)
	cai  		(mM)

	caix = 3
									





        Ra = 0.2
        Rb = 0.2
        Rc = 0.5
        Rd = 0.5

	dt		(ms)
	celsius		(degC)
	temp = 23	(degC)		
	q10  =  3			

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
	ntau 		(ms)
        hinf
        htau	
	tadj
}
 

STATE { n h }

INITIAL { 
	rates(cai)
	n = ninf
        h = hinf
}

BREAKPOINT {
        SOLVE states
	gk = tadj*gbar*(n^(0.2))*h
        
	
        ik = gk * (v - ek)
} 

LOCAL nexp, hexp

PROCEDURE states() {    
        rates(cai)      
        n = n + nexp*(ninf-n)
        h = h + hexp*(hinf-h)

        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(cai(mM)) {  

        LOCAL tinc

        a = Ra * cai^caix
        b = Rb
        ntau = 1/(a+b)
      	ninf = a*ntau
        
        a = Rc * cai^(-5)
        b = Rd
        htau = 1/(a+b)
	hinf = a*htau

        tadj = q10^((celsius - temp)/10)

        tinc = -dt * tadj
        nexp = 1 - exp(tinc/ntau)
	hexp = 1 - exp(tinc/htau)	
}