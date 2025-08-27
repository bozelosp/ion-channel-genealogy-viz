NEURON {
	SUFFIX capq
	USEION ca READ eca WRITE ica
	RANGE gca, gbar, i
	RANGE minf, mtau 
	GLOBAL midv, mslope, k1, k2, v1, v2, temp
}


UNITS { 
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	
	eca = 110 (mV) 
	gbar = 1 (S/cm2) 	
	midv = -3.9 (mV)
	mslope = 7.1 (mV)
	k1= 0.22 
	k2= 0.14
	v1= 31.5 (mV)
	v2= 8.6 (mV)
	delay=0.2 (ms)
	temp=0
}

ASSIGNED {
	
	v (mV)
	i 	(mA/cm2)
	ica 	(mA/cm2)
	gca	(S/cm2)
	minf
	mtau (ms)
}

STATE { m }

INITIAL { 
	rates(v)
	m = minf
}

BREAKPOINT {
		SOLVE states METHOD cnexp
        	gca = gbar*m
		i = gca*(v - eca) 
		ica = i
}

DERIVATIVE states {  
		rates(v)
        m' = (minf - m)/mtau
}

PROCEDURE rates(v (mV)) 
{
	minf = 1/ ( 1 + exp( (midv-v)/mslope ) )
	mtau = ((2.3^-temp)/(k1*exp(v/v1)+k2*exp(-v/v2)))+delay
}