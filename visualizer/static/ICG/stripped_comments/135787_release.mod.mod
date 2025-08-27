INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX rel
	USEION ca READ cai WRITE cai
	RANGE T,FA,CA,Fmax,Ves,b,u,k1,k2,k3,nt,kh
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {

	Ves = 0.1 	(mM)		
	Fmax = 0.001	(mM)		
	b = 1e16 	(/mM4-ms)	
	u = 0.1  	(/ms)		
	k1 = 1000   	(/mM-ms)	
	k2 = 0.1	(/ms)		
	k3 = 4   	(/ms)		
	nt = 10000			
	kh = 10  	(/ms)		
}

ASSIGNED {
}

STATE {
	FA	(mM)
	VA	(mM)
	T	(mM)
	cai	(mM) 
}

INITIAL {
	FA = 0
	VA = 0
	T = 0
	cai = 1e-8
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit 
}

LOCAL bfc , kfv

DERIVATIVE state {

	bfc = b * (Fmax-FA-VA) * cai^4
	kfv = k1 * FA * Ves

	cai'	= - bfc + 4 * u * FA
	FA'	= bfc - u * FA - kfv + k2 * VA
	VA'	= kfv - (k2+k3) * VA
	T'	= nt * k3 * VA - kh * T
}