UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	gkdrbar=.003 (mho/cm2)
	vhalfn=13	(mV)
	a0n=0.02	(/ms)
	zetan=-3	(degC/mV)
	gmn=0.7		(1)
	nmax=2		(ms)
	q10=1		(1)
}


NEURON {
    THREADSAFE

	SUFFIX kdr
	USEION k READ ek WRITE ik
	RANGE gkdr,gkdrbar

}

STATE {
	n
}

ASSIGNED {
	v		(mV)
	ek		(mV)
	celsius	(degC)
	ik		(mA/cm2)
	ninf	(1)
	gkdr	(mho/cm2)
	taun	(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gkdrbar*n
	ik = gkdr*(v-ek)

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16(degC)+celsius))) 
}

DERIVATIVE states {     
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { 
        LOCAL a,qt
        qt=q10^((celsius-24(degC))/10(degC))
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
	if (taun<nmax) {taun=nmax}
}