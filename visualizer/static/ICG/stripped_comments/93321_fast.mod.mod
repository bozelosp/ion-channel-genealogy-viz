NEURON {
	SUFFIX F
	USEION ca READ ica
	RANGE F, M, H
}

UNITS {
	(mA) = (milliamp)
}

PARAMETER {
        G = 10           
	tau_M = 0.5 (ms) 
	tau_H = 1.5 (ms) 
	Z_M = 14.2       
	Z_H = 9.8        
}

ASSIGNED {
	Mbar (1)  
	Hbar (1)  
	ica (mA/cm2)
	F (1)
}

STATE {	M H }

BREAKPOINT {
	SOLVE states METHOD cnexp
	F = G * M * M * H    
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