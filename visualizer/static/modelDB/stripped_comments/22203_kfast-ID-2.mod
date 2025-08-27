INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kfast
	USEION k READ ek WRITE ik
	RANGE  a, b, gkfast, gbar
	RANGE  ainf, taua, binf, taub

	GLOBAL a0, a1, a2, a3, a4, a5, a6
	GLOBAL b0, b1, b3, b4, b5, b6
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0   	(pS/um2)	
	vshift = 0	(mV)		
									
	a0   =  3.4  	(1/ms)		
	a1   =  10 	(mV)		
	a2   =  30 	(mV)		

	a3   = 0.7 	(1/ms)		
	a4   = 6    	(mV)		
	a5   = 180  	(mV)		
	a6   = -0.45 	(1/ms)		

	b0   = 0.0001          (1/ms)		
	b1   = 0.06	(1/mV)		
	


	b3   = 0.12    (1/ms)		
	b4   = -50      (mV)
	b5   = 8         (mV)
	b6   = 0.002  (1/ms)
	

	temp = 21	(degC)		
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
	gkfast		(pS/um2)
	ek		(mV)
	ainf 		
	binf
	taua (ms)	
	taub (ms)	
	tadj
}
 

STATE { a b }

INITIAL { 
	trates(v+vshift)
	a = ainf
	b = binf 
}

BREAKPOINT {
        SOLVE states
        gkfast = tadj*gbar*a*a*a*a*b
	ik = (1e-4) * gkfast * (v - ek)
} 

LOCAL aexp, bexp, z 

PROCEDURE states() {   		
        trates(v+vshift) 	
        a = a + aexp*(ainf-a)
        b = b + bexp*(binf-b)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(v) { 
                      
    LOCAL tinc

    TABLE ainf, binf, aexp, bexp, taua, taub
	DEPEND dt, celsius, temp, q10, a0, a1, a2, a3, a4, a5,a6, b0, b1, b3, b4, b5, b6
	FROM vmin TO vmax WITH 199

	rates(v)

        tadj = 1 
        tinc = -dt * tadj

        aexp = 1 - exp(tinc/taua)
        bexp = 1 - exp(tinc/taub)
}



PROCEDURE rates(vm) {  

	LOCAL alphaA, betaA,AlphaI,BetaI
	
	alphaA=a0/(1+exp(-(vm-a1)/a2))
	betaA=a3*exp(-(vm-a4)/a5)+a6
    	taua = 1/(alphaA+betaA)
	ainf = alphaA/(alphaA+betaA)


	AlphaI=b0*exp(-b1*vm)
	BetaI=b3/(1+exp(-(vm-b4)/b5))+b6
    	taub = 1/(AlphaI+BetaI)
	binf = AlphaI/(AlphaI+BetaI)
}