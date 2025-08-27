NEURON { 
	SUFFIX TNC 
	NONSPECIFIC_CURRENT i
	RANGE gbar, i, eTNC
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
    gbar = 1e-5 (siemens/cm2)
} 

ASSIGNED {
	v (mV)
	eTNC (mV) 
    i (mA/cm2)
} 
 
BREAKPOINT { 
	i = gbar * (v - eTNC)
}