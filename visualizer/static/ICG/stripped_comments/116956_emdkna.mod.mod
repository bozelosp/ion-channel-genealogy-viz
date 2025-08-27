NEURON {
	SUFFIX emdkna
	USEION k READ ek WRITE ik
	USEION na READ ina
	RANGE gk, gbar, i
	RANGE ninf, ntau 
	GLOBAL slope,taumax
}

UNITS { 
  (mV) = (millivolt)
  (mA) = (milliamp)
  (uA) = (microamp)
  (S) = (siemens)
}

PARAMETER {
	ek = -20 (mV) 
	gbar = 0.002 (S/cm2) 	
	
	slope = 7 (uA/cm2)
	taumax = 3 (ms)
	ntau = 3 (ms) 
}

ASSIGNED {
	v (mV) 
	i 	(mA/cm2)
	ik 	(mA/cm2)
	gk		(S/cm2)
	ninf
	ina (mA/cm2) 
}

STATE { n }

INITIAL { 
	rates(ina)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n 
	i = gk * (v - ek)
	ik = i
}

DERIVATIVE states {  
	rates(ina)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(ina (mA/cm2)) {
	ninf = 1 / ( 1 + exp( -ina/((1e-3)*slope) ) )  
}