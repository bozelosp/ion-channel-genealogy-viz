NEURON { 
	SUFFIX GrC_Lkg2 
	NONSPECIFIC_CURRENT il
	RANGE egaba, ggaba , i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
		ggaba = 2.17e-5 (mho/cm2)
	      egaba = -65 (mV)
} 

ASSIGNED { 
       v (mV) 
       il (mA/cm2) 
	 i (mA/cm2) 
}

BREAKPOINT { 
	il = ggaba*(v - egaba) 
	i =il
}