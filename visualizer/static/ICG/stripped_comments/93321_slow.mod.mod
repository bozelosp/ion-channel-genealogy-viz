NEURON {
	SUFFIX S
	USEION ca READ ica
	RANGE S
}

UNITS {
	(mA)	=	(milliamp)
}

PARAMETER {
	G = 3           
	tau_M = 50 (ms) 
	tau_H = 60 (ms) 
	Z_M = 7.2       
	Z_H = 2.8       
}

ASSIGNED {
	Mbar (1)  
	Hbar (1)  
	ica (mA/cm2)
	S (1)
}

STATE {	M H }

BREAKPOINT {
	SOLVE states METHOD cnexp
	S = G * M * M * H    
}

INITIAL {
	
	rates(ica)
	M = Mbar
	H = Hbar
}
DERIVATIVE states {
	rates(ica)
	M' = (Mbar - M)/tau_M
	H' = (Hbar - H)/tau_H
}

PROCEDURE rates(ica (mA/cm2)) {
	UNITSOFF
	Mbar = 1/(1+exp( Z_M + 1e3*ica ))  
	Hbar = 1/(1+exp( -Z_H - 1e3*ica )) 
	UNITSON
}