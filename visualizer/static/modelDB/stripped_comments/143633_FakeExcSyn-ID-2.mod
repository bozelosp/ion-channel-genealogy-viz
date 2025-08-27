INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS FakeExcSyn
	POINTER pre
	RANGE lastrelease, dummy, lastspike, spike, delay, amp
	NONSPECIFIC_CURRENT i
	GLOBAL dur, Prethresh, Deadtime
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	dur	= 0.1	(ms)		
	amp = 2  (nA) 
	Prethresh = -20 			
	Deadtime = 1	(ms)		
	delay = 0 (ms)			
							
							
}


ASSIGNED {
	i 		(nA)		
	pre 				
	lastrelease	(ms)		
	lastspike (ms)		
	dummy
	spike
}

INITIAL {
	lastrelease = -9e9
	lastspike = 0
	dummy = 0
	spike = 0
}

BREAKPOINT {
	SOLVE release
	i = -dummy*amp
}

PROCEDURE release() { LOCAL q, r
	

	q = ((t - lastrelease) - dur)		
	r = (t - lastspike)					
		
						
	if (q > Deadtime && spike == 0) {
		
		if (pre > Prethresh) {		
			spike = 1
			lastspike = t
		}
				
	} else if (q > Deadtime && spike == 1) {
		
		if (t >= lastspike + delay) {
			dummy = 1			
			lastrelease = t
		}
		
	} else if (q < 0) {			
	
		
	
	} else if (dummy == 1) {			
		dummy = 0
		spike = 0
	}
}