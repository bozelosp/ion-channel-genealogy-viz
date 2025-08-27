INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX rel
	USEION ca READ ica, cai WRITE cai
	RANGE T,FA,CA,Fmax,Ves,b,u,k1,k2,k3,nt,kh

	RANGE depth,kt,kd,cainf,taur
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)

	(molar) = (1/liter)			

	(um)	= (micron)

	(msM)	= (ms mM)

}


CONSTANT {
	FARADAY = 96489		(coul)		

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

	depth	= .1	(um)		
	taur	= 700	(ms)		
	cainf	= 1e-8	(mM)
	kt	= 1	(mM/ms)		
	kd	= 5e-4	(mM)		
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
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

	cai = kd
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

LOCAL bfc , kfv

DERIVATIVE state {

	bfc = b * (Fmax-FA-VA) * cai^4
	kfv = k1 * FA * Ves


	FA'	= bfc - u * FA - kfv + k2 * VA
	VA'	= kfv - (k2+k3) * VA
	T'	= nt * k3 * VA - kh * T


	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	


	drive_pump = -kt * cai / (cai + kd )		



	cai'= -bfc+4*u*FA + drive_channel + drive_pump + (cainf-cai)/taur

}