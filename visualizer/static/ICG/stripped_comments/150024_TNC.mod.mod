NEURON { 
	SUFFIX TNC 
	NONSPECIFIC_CURRENT i
	RANGE gbar, i
	GLOBAL eh
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
    gbar = 1e-5 (siemens/cm2)
} 

ASSIGNED {
        eh (mV)
	v (mV)
	
    i (mA/cm2)
} 
 
BREAKPOINT { 
	i = gbar * (v - eh)
}