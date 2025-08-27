NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE n, gbar,ik
	GLOBAL Ra, Rb, ninf, ntau
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)


} 

PARAMETER {
	gbar = 0.03     (mho/cm2)

	tha  = -30	(mV)		
	qa   = 9	(mV)		
	Ra   = 0.001	(/ms)		
	Rb   = 0.001	(/ms)		
	temp = 23	(degC)		
	q10  = 2.3			
	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	celsius		(degC)
	v 		(mV)
	ik 		(mA/cm2)
	ek		(mV)
	ninf
	ntau 		(ms)
	tadj
}
 

STATE { 
	n 
}

INITIAL {
        tadj = q10^((celsius - temp)/10(degC))  
	rates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp

	ik = (1e-4) *tadj* gbar*n * (v - ek)
}


DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/ntau
}


PROCEDURE rates(v(mV)) {  
                      
	LOCAL a,b
        a = Ra * (v - tha)*1(/mV) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha)*1(/mV) / (1 - exp((v - tha)/qa))
        ntau = 1/(a+b)/tadj
	ninf = a*ntau
}