NEURON { 
	SUFFIX GRC_LKG1 
	NONSPECIFIC_CURRENT il
	RANGE el, gl,i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gl = 5.68e-5 (mho/cm2)
	celsius = 30 (degC)
	el =  -16.5 (mV) 
	
} 

ASSIGNED { 
	il (mA/cm2) 
	i (mA/cm2) 
}
  
BREAKPOINT { 
	il = gl*(v - el) 
	i = il
}