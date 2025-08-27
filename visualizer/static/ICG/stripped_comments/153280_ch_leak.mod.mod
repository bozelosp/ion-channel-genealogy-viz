VERBATIM
#include <stdlib.h> 
ENDVERBATIM

UNITS {
	(mA) =(milliamp)
	(mV) =(millivolt)
}
 
NEURON { 
	SUFFIX ch_leak 
	NONSPECIFIC_CURRENT i
	RANGE gmax, e, i
	RANGE myi
    THREADSAFE
}
 
PARAMETER {
	g (mho/cm2)		
	gmax (mho/cm2)		
	e (mV)			
}

ASSIGNED {	
	v (mV) 			
					
					
	i (mA/cm2)		
	myi (mA/cm2)
} 

BREAKPOINT {
	g = gmax
	i = g*(v-e)	
	myi = i
}