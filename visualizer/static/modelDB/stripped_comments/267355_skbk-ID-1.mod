UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX bks
	USEION lca READ lcai
USEION nca READ ncai
USEION tca READ tcai
	USEION k WRITE ik
	RANGE gkbar,gk,zinf,ik
	GLOBAL bkcoef,bkexp,cahco, carco,tin, tinsh, calco
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.08	(mho/cm2)	
	 ncai 	(mM)
	 lcai (mM)
	 tcai (mM)

	ek  = -85	(mV)
	dt		(ms)
	tin=0.02
	bkcoef=0.023
	bkexp=7.5
	cahco=1
	carco=0
	calco=1
	tinsh=400
}


ASSIGNED {
	ik		(mA/cm2)
	minf
	mexp
	zinf
	zexp
	gk
}

STATE {	m z }		

BREAKPOINT {
	SOLVE state

	ik = gkbar*1000*m*z*z*(v - ek)
}






PROCEDURE state() {	
	rate(v, lcai,ncai,tcai)
	m = m + mexp*(minf - m)
	z = z + zexp*(zinf - z)
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	rate(v, lcai,ncai, tcai)
	m = minf
	z = zinf
}

FUNCTION alp(v (mV), lcai (mM), ncai (mM), tcai(mM)) (1/ms) { 
	alp = 0.4/((lcai*calco+ncai*cahco+tcai*carco))
}


FUNCTION bet(v (mV)) (1/ms) { 
	bet = 0.11/exp((v-55)/14.9)
}

PROCEDURE rate(v (mV), lcai (mM),ncai (mM), tcai (mM)) { 
	LOCAL a,b,tinca
	a = alp(v,lcai,ncai, tcai)
	zinf = 1/(1+a)
	
	zexp = (1 - exp(-dt/4))
	b = bet(v)
	minf = 8.5/(7.5+b)
	mexp = (1 - exp(-dt*(bkexp+b)))
}