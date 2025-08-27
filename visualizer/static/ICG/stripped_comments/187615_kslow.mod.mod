INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kslow
	USEION k READ ek WRITE ik
	RANGE  a, b, b1,gkslow, gbar, ik
	RANGE  ainf, taua, binf, taub,taub1
	GLOBAL a0, a1, a2, offma, sloma, offmb, slomb
	GLOBAL b0, b11, b2, offht1, sloht2, offht2
	GLOBAL bb0,bb1,bb2, offht3, offht4, sloht3, sloht4
	GLOBAL offh, sloh
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 1   	(pS/um2)	
	vshift = 0	(mV)		
								
	offh = -58	(mV)		
	sloh   = 11  (mV)		
		
	a0   = 192.3076923 (ms mV)		
	a1   = 51.5995872 (ms)
	a2   = 188.6792453 (ms)	
	offma   = 11.1	(mV)
	sloma   = 13.1	(mV)
	offmb   = -1.27	(mV)
	slomb   = 71    (mV)
	
	b0   = 360	(ms)			
	b11   = 1010	(ms)	
	b2   = 23.7     (ms/mV)
	offht1   = -54      (mV)
	offht2   = -75	(mV)	
	sloht2   = 48	(mV)	

	bb0 = 2350	(ms)			
	bb1 = 1380	(ms)
	bb2 = -210  (ms)
	offht3 = 0	(mV)
	offht4 = 0	(mV)
	sloht3 = 89.4454383 (mV)
	sloht4 = 32.67973856 (mV)

	temp = 21	(degC)		
	q10  = 2.3			

	v 		(mV)
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
	gkslow		(pS/um2)
	ek		(mV)
	ainf 		
	binf
	taua (ms)	
	taub (ms)
	taub1 (ms)	
	tadj
}
 

STATE {a b b1}

INITIAL { 
	rates(v+vshift)
	a = ainf
	b = binf 
	b1= binf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gkslow = gbar*a*a*(0.5*b+0.5*b1)
	  ik = (1e-4) * gkslow * (v - ek)
} 

LOCAL aexp, bexp,b1exp, z 

DERIVATIVE states {   		
        rates(v+vshift) 	
        a'  = (ainf-a)/taua
        b'  = (binf-b)/taub
	b1' = (binf-b1)/taub1
}


PROCEDURE rates(vm) {  

	LOCAL alpha, beta
	
	tadj = q10^((celsius - temp)/10)
	
	alpha=tadj/a0*(vm-offma)/(1-exp(-(vm-offma)/sloma))
	beta=tadj/a1*exp(-(vm-offmb)/slomb)-1/a2

	taua=1/(alpha+beta)
	ainf = alpha/(alpha+beta)
	
	taub = b0 + (b11+b2*(vm-offht1))*exp(-(vm-offht2)*(vm-offht2)/(sloht2*sloht2))
    	taub1=bb0+bb1*exp(-(vm-offht3)/sloht3)+bb2*exp(-(vm-offht4)/sloht4)
	binf = 1/(1+exp((vm-offh)/sloh))
}