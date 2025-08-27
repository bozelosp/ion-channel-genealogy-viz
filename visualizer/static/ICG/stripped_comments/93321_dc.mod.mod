NEURON {
	SUFFIX D
	USEION ca READ ica
	RANGE D
}

UNITS {
	(mA) = (milliamp)
}

PARAMETER {
	G = 1            
	tau_M = 500 (ms) 
	Z_M = 3          
}

ASSIGNED {
	Mbar (1)  
	ica (mA/cm2)
	D (1)
}

STATE {	M (1)}

BREAKPOINT {
	SOLVE states METHOD cnexp
	D = G * M * M    
}

INITIAL {
	
	rates(ica)
	M = Mbar
}
DERIVATIVE states {
	rates(ica)
	M' = (Mbar - M)/tau_M
}

PROCEDURE rates(ica (mA/cm2)) {
	UNITSOFF
	Mbar = 1/(1+exp( Z_M + 1e3*ica )) 
	                              
	UNITSON
}