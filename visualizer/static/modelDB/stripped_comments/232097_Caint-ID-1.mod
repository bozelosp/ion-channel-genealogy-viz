NEURON {
	SUFFIX Cacon
	USEION ca READ ica, cai WRITE cai
	GLOBAL tauca, A, camin
}

UNITS { 
 (mA) = (milliamp) 
 (mV) = (millivolt) 
 (molar) = (1/liter) 
 (mM) = (millimolar)
 (uM) = (micromolar) 
} 

PARAMETER {  
 dt (ms) 
 tauca = 800 (ms) 

 A = 0.2  
 camin = 1e-8 (mM)  

} 

STATE {
	cai		(mM) 
}

INITIAL {
	cai = camin
}

ASSIGNED { 
 ica (mA/cm2)    
} 

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	 cai' = -A*ica - (cai-camin)/tauca 
}