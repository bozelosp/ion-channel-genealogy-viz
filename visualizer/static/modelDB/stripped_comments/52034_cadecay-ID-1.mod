INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}

NEURON {
	SUFFIX cadecay
	USEION ca READ ica WRITE cai
	RANGE taucaremov
	GLOBAL cainit
}


UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)

	FDY	= 9.6520	( 10000 coulomb)
}


PARAMETER {
	taucaremov	= 20	(ms)
	diam 		= 1	(um)
	ica			(mA/cm2)
	cainit		= 5e-5	(mM)
}

STATE {
	cai		(mM) 
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	cai' = ((cainit - cai)/taucaremov) + (-ica * 4/(diam*FDY))


}

INITIAL {
    cai = cainit
}