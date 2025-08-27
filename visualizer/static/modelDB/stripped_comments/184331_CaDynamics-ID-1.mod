NEURON	{
	SUFFIX CaDynamics
	USEION ca READ ica WRITE cai
	RANGE decay, gamma, minCai, depth
}

UNITS	{
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um)	= (micron)
}

PARAMETER	{
	gamma = 0.05 
	decay = 80 (ms) 
	depth = 0.1 (um) 
	minCai = 1e-4 (mM)
}

ASSIGNED	{ica (mA/cm2)}

INITIAL {
	cai = minCai
}

STATE	{
	cai (mM)
}

BREAKPOINT	{ SOLVE states METHOD cnexp }

DERIVATIVE states	{
	cai' = -(10000)*(ica*gamma/(2*FARADAY*depth)) - (cai - minCai)/decay
}