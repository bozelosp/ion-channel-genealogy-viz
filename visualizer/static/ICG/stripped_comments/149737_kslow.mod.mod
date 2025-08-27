INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kslow
	USEION k READ ek WRITE ik
	RANGE  a, b, b1,gkslow, gbar
	RANGE  ainf, taua, binf, taub,taub1
	GLOBAL a0, a1, a2, a3, a4, a5, a6
	GLOBAL b0, b11, b2, b3, b4, b5
	GLOBAL bb0,bb1,bb2,bb3,bb4
	GLOBAL v05a, za, v05b, zb
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 1.0   	(S/cm2)	
	vshift = 0	(mV)		
								
	v05a = -14.3	(mV)		
	za   =  14.6	(mV)		
	v05b = -58	(mV)		
	zb   = -11  (mV)		
		
	a0   =  0.0052  (1/ms 1/mV)		
	a1   = 11.1 	(mV)			
	a2   = 13.1	(mV)				
	a3   = 0.02    (1/ms)		
	a4   = -1.27	(mV)			
	a5   = 71    (mV)
	a6   = -0.005 (1/ms)	
	
	b0   = 360	(ms)			
	b11   = 1010	(ms)		
	b2   = -75	(mV)			
	b3   = 48	(mV)			
	b4   = 23.7     (ms/mV)
	b5   = -54      (mV)

	bb0 = 2350	(ms)			
	bb1 = 1380	(ms)
	bb2 = 0.011 (mV)
	bb3 = -210  (ms)
	bb4 = 0.03 (mV)

	temp = 21	(degC)		
	q10  = 2.3			

	ek   (mV)
	ekkai  = -80 (mV)
	v 		(mV)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gkslow		(S/cm2)
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
	  ik = gkslow * (v - ekkai)
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
	
	alpha=a0*(vm-a1)/(1-exp(-(vm-a1)/a2))
	beta=a3*exp(-(vm-a4)/a5)+a6

	taua=1/(alpha+beta)
	ainf = alpha/(alpha+beta)
	
	taub = b0 + (b11+b4*(vm-b5))*exp(-(vm-b2)*(vm-b2)/(b3*b3))
    	taub1=bb0+bb1*exp(-bb2*vm)+bb3*exp(-bb4*vm)
	binf = 1/(1+exp(-(vm-v05b)/zb))
}