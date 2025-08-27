TITLE Calcium dependent potassium channel (small conductance SK) from callaway (2000)
 

UNITS {
        (molar) = (1/liter)
        (pA) =  (picoamp)
	(mV) =	(millivolt)
        (pS)  =  (picosiemens)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX kca
	USEION ca READ cai
	USEION k WRITE ik
	RANGE  gbar,k_half, oinf, ik, g, n
 
}


PARAMETER {
        dt  (ms)
        cai (mM)
        celsius = 32  (degC)
        gbar =0.1 (pS/microm2)
        ek = -90.0      (mV)
        k_half = 0.00019   (mM)
    	 n=4.0 (1)
        
        
}

ASSIGNED { 
           ik		(mA/cm2)
           oinf           
		    g	(pS/microm2)
}


BREAKPOINT {
        oinf =1/(1+pow(k_half/cai,n))
        		g=gbar*oinf
	ik = (0.0001)*g*(v - ek)
}

