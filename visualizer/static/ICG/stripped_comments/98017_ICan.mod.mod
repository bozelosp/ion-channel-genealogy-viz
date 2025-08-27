UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(mM) = (milli/liter)
}

NEURON {
        SUFFIX ICan
	USEION ca READ cai
        NONSPECIFIC_CURRENT i
        RANGE gbar, g, e, minf, sum
}

PARAMETER {
        gbar 	= 0.000025	(mho/cm2)	<0,1e9>
	e   	= -20 		(mV)
	Beta   	= 0.0004	(/ms)
	Alpha	= 5000		(/ms-mM2)
}
 
ASSIGNED {
        v 	(mV)
	g 	(mho/cm2)
        i 	(mA/cm2)
	minf	(1)
	taum	(ms)
	cai	(mM)
	sum	(/ms)
}

STATE { m }

INITIAL{ 
	m = minf
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
	rates(cai)
        g = gbar*m*m
	i = g * ( v - e )
}

DERIVATIVE state {
	m' = ( minf - m ) / taum
}

PROCEDURE rates( cai (mM) ) { 
	sum = Afac(cai) + Beta
	taum = 1 / sum
	minf = Afac(cai) / sum
}

FUNCTION Afac( cai (mM) ) (/ms) {
	Afac = Alpha * cai * cai
}