UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
	(S) = (siemens)
}

NEURON {
	SUFFIX kaslow
	USEION k READ ek WRITE ik
	RANGE gbar, gka, ik
	RANGE ninf, ntau, hinf, h1tau, h2tau
	THREADSAFE
}
 
PARAMETER {
        gbar = 0.00136 (S/cm2)
		
		A1 = 25.46
		B1 = 67.41
		xc1 = 50
		yc1 = 21.95
		
		A2 = 200
		B2 = 587.4
		xc2 = 0
		yc2 = 47.77
		
		AN =  1.1972
		BN = 2.56
		xcN = 60
		ycN = 45.75992
		
		vh1 = -40.8 (mV)
		vh2 = 74.2 (mV)
		
		frac1 = 0.3
		frac2 = 0.7
}
 
STATE {
        n 	
		h1  
		h2	
}
 
ASSIGNED {
	v (mV)
	ek (mV)
	gka (S/cm2)
    ik (mA/cm2)
	    
    ninf
	ntau (ms)
	hinf
	h1tau (ms)
	h2tau (ms)
}
 
BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gbar*n*(h1*frac1+h2*frac2)
	ik = gka*(v - ek)
}
 
INITIAL {
	rates(v)
	
	n = ninf
	h1 = hinf
	h2 = hinf
}

DERIVATIVE states {  
	rates(v)
	
	n'  = (ninf-n)/ntau
	h1' = (hinf-h1)/h1tau
	h2' = (hinf-h2)/h2tau
}
 
PROCEDURE rates(v(mV)) { 
        UNITSOFF   
		ninf = (1/(1+exp((vh1-v)/9.5)))  		
		ntau = AN + BN*exp(-2*((v+xcN)/ycN)^2)	
		
		hinf  = 1/(1 + exp((v+vh2)/9.6))  		
		h1tau = A1 + B1*exp(-2*((v+xc1)/yc1)^2) 
		h2tau = A2 + B2*exp(-((v-xc2)/yc2)^2)	
}
 
UNITSON