NEURON {
	SUFFIX emdk
	USEION k READ ek WRITE ik
	RANGE gk, gbar, i
	RANGE ninf, ntau 
	GLOBAL nmidv, nslope, ntaumax, nmidvdn, nslopedn, nmidvup, nslopeup
}


UNITS { 
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	
	ek = -20 (mV) 
	gbar = 0.001 (S/cm2) 	
	nmidv = 14 (mV)
	nslope = 11 (mV)
	ntaumax = 50.2 (ms)
	nmidvdn = 25 (mV)
	nslopedn = -30 (mV)
	nmidvup = 28 (mV)
	nslopeup = 32 (mV)
}

ASSIGNED {
	
	v (mV)
	i 	(mA/cm2)
	ik 	(mA/cm2)
	gk	(S/cm2)
	ninf
	ntau (ms)
}

STATE { n }

INITIAL { 
	rates(v)
	n = ninf
}

BREAKPOINT {
		SOLVE states METHOD cnexp
        gk = gbar*n*n*n*n 
		i = gk * (v - ek) 
		ik = i
}

DERIVATIVE states {  
		rates(v)
        n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) 
{
	ninf = 1/ ( 1 + exp( (nmidv-v)/nslope ) )
	ntau = ntaumax / ( exp( (nmidvdn-v)/nslopedn ) + exp( (nmidvup-v)/nslopeup ) )
	
}