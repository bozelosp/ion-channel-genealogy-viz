INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON 
{
	POINT_PROCESS AMPA
	RANGE C, g, gmax, lastrelease, TRise, tau
	NONSPECIFIC_CURRENT i
	RANGE Cmax, Cdur, Alpha, Beta, Erev, Deadtime
}

UNITS 
{
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER 
{
	TRise  	= 2 (ms)
	tau    	= 2(ms)
	Cmax	= 1	(mM)		
	Erev	= 0	(mV)		
	Deadtime = 1	(ms)		
	gmax	= 0		(umho)		
}


ASSIGNED 
{
	Alpha	(/ms mM)	
	Beta	(/ms)		
	Cdur	(ms)		
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C		(mM)		
	lastrelease	(ms)		
}

STATE
{
	R				
}

INITIAL 
{
	R = 0
	C = 0
	lastrelease = -1000
	Cdur=TRise
	Beta=1/tau
	Alpha=1/Cdur - Beta
}

BREAKPOINT 
{
	SOLVE states METHOD cnexp
	g = (gmax * R * (Alpha+Beta)) / (Alpha*(1-1/exp(1)))
	i = g*(v - Erev)
}

DERIVATIVE states
{
	evaluateC() 	
	R'=Alpha * C * (1-R) - Beta * R
}

PROCEDURE evaluateC()
{
	LOCAL q
	q = ((t - lastrelease) - Cdur)		
	if (q >= 0 && q <= Deadtime && C == Cmax) {	
		C = 0.
	}
}

NET_RECEIVE (weight (umho)) 
{ 
	LOCAL q
	q = ((t - lastrelease) - Cdur)		



	if (q > Deadtime) {
		C = Cmax			
		lastrelease = t
	} 
}