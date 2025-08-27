NEURON {
	SUFFIX gh
	
	
	NONSPECIFIC_CURRENT i
	GLOBAL eh
	RANGE ghbar, ik, ina,htau, half, slp,inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	ghbar =.001 (mho/cm2) <0,1e9>
	htau = 50 (ms)
	half=-80 (mV)
	slp=8 (mV)
	
	
}
STATE {
	n
}
ASSIGNED {
	
	
	i (mA/cm2)
	eh (mV)
	inf
}

INITIAL {
	rate(v)
	n = inf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	
	
	i = ghbar*n*(v-eh)
}

DERIVATIVE states {	
	rate(v)
	n' = (inf - n)/htau
}
UNITSOFF

PROCEDURE rate(v(mV)) {	
	TABLE inf DEPEND half,slp FROM -100 TO 100 WITH 200
		inf = 1/(1+exp((v-half)/slp))
}
UNITSON