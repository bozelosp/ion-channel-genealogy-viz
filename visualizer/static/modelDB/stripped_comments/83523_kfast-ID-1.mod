INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kfast
	USEION k READ ek WRITE ik
	RANGE  a, b, gbar
	RANGE  ainf, taua, binf, taub,gkfast
	GLOBAL a0, a1, a2, a3, a4, a5, a6
	GLOBAL b0, b1, b2, b3, b4, b5
	GLOBAL v05a, za, v05b, zb
	GLOBAL q10, temp, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0   	(pS/um2)	
	vshift = 0	(mV)		
								
	v05a = -56	(mV)		
	za   =  26	(mV)		
	v05b = -75	(mV)		
	zb   = -11  (mV)		
		
	a0   =  0.096  	(1/ms)		
	a1   =  0.18 	(1/ms)		
	a2   =  0.02 	(1/mV)		

	a3   = 0.104 	(1/ms)		
	a4   = 1.64    	(1/ms)		
	a5   = -16.88  	(mV)			
	a6   = 9 		(mV)			
	
	b0   = 956	(ms)			
	b1   = 16.6	(ms/mV)	
	b2   = 0.076	(ms/mV^2)	

	b3   = 6.5	(ms)			
	b4   = 0.00082     (ms)
	b5   = -0.147      (1/mV)


	temp = 33	(degC)		
	q10  = 2.3				

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	ek		(mV)
	gkfast	(pS/um2)
	ainf 		
	binf
	taua (ms)	
	taub (ms)
}
 

STATE { a b }

INITIAL { 
	rates(v+vshift)
	a = ainf
	b = binf 
}

BREAKPOINT {
      SOLVE states METHOD cnexp
	gkfast = gbar*a^4*b
	ik = (1e-4) * gkfast * (v - ek)
} 


DERIVATIVE states {   		
        rates(v+vshift) 	
        a' = (ainf-a)/taua
        b' = (binf-b)/taub
}


PROCEDURE rates(vm  (mV)) {  

	LOCAL alphaA, betaA
	TABLE ainf, binf, taua, taub
	DEPEND dt, celsius, temp, q10, a0, a1, a2, a3, a4, a5, a6, b0, b1, b2, b3, b4,b5, v05a, za, v05b, zb
	FROM vmin TO vmax WITH 199

	ainf = 1/(1+exp(-(vm-v05a)/za))
	alphaA=a0+a1*exp(-a2*vm)
	betaA=a3+a4/(1+exp(-(vm-a5)/a6))
    	taua = betaA*(1-heav(vm+40))+alphaA*heav(vm+40)

	
	taub = (b0+b1*vm+b2*vm^2)*(1-heav(vm+80))+(b3+b4*exp(b5*vm))*heav(vm+80)
	binf = 1/(1+exp(-(vm-v05b)/zb))
}



FUNCTION heav(y) {
	LOCAL z
	IF (y<=0) {     
	z=0 
		}
		ELSE {
		z=1
		}
		heav=z 
}