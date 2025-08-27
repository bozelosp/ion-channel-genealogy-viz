NEURON {
	SUFFIX emdna
	USEION na READ ena WRITE ina 
	RANGE gna, gbar, i
	RANGE minf, hinf, mtau, htau 
	GLOBAL hmidv, hslope, htaumax, hmidvdn, hslopedn, hmidvup, hslopeup
	GLOBAL mmidv, mslope, mtaumax, mmidvdn, mslopedn, mmidvup, mslopeup
}

UNITS {
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	ena = 100 (mV) 
	gbar = 0.003 (S/cm2) 
	
	mmidv = -1 (mV)
	mslope = 6 (mV)
	mtaumax = 1.8 (ms)
	mmidvdn = -31 (mV)
	mslopedn = -25 (mV)
	mmidvup = -60 (mV)
	mslopeup = 11 (mV)
	
	hmidv = -11 (mV)
	hslope = -8 (mV)
	htaumax = 9 (ms)
	hmidvdn = -4 (mV)
	hslopedn = -12 (mV)
	hmidvup = -16 (mV)
	hslopeup = 14 (mV)
}

ASSIGNED {
	v (mV)
	i (mA/cm2)
	ina (mA/cm2)
	gna	(mho/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
}

STATE { m h }

INITIAL { 
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
		SOLVE states METHOD cnexp
        gna = gbar*m*m*m*h
		i = gna * (v - ena)
		ina = i
} 

DERIVATIVE states {  
		rates(v)
        m' = (minf - m)/mtau
        h' = (hinf - h)/htau
}

PROCEDURE rates(v (mV)) {
	
	minf = 1/ ( 1 + exp( (mmidv-v)/mslope ) )
	mtau = mtaumax / ( exp( (mmidvdn-v)/mslopedn ) + exp( (mmidvup-v)/mslopeup ) )
	
	
	
	hinf = 1/ ( 1 + exp( (hmidv-v)/hslope ) )
	htau = htaumax / ( exp( (hmidvdn-v)/hslopedn ) + exp( (hmidvup-v)/hslopeup ) )
	
}