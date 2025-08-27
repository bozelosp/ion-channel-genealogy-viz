UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt (ms)
	v (mV)
        ek = -90	(mV)	
	celsius = 24	(degC)
	gkdrbar=.003 (mho/cm2)
        ikmax = 0.3 (mA/cm2)
        vhalfn=13   (mV)
        a0n=0.02      (/ms)
        zetan=-3    (1)
        gmn=0.7  (1)
	nmax=2  (ms)
	q10=1
        nscale=1
}


NEURON {
	SUFFIX kdr
	USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar,ik
	RANGE ninf,taun
        GLOBAL nscale
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkdr (mho/cm2)
        taun (ms)
}

INITIAL {
        rates(v)
        n=ninf
	gkdr = gkdrbar*n
	ik = gkdr*(v-ek)
}        

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gkdrbar*n
	ik = gkdr*(v-ek)

}

DERIVATIVE states {
        rates(v)
        n' = (ninf-n)/taun
}

FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4(degC/mV)/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4(degC/mV)/(8.315*(273.16+celsius))) 
}

LOCAL facn










PROCEDURE rates(v (mV)) { 
        LOCAL a,qt
        qt=q10^((celsius-24)/10(degC))
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
	if (taun<nmax) {taun=nmax}
        taun=taun/nscale
        facn = (1 - exp(-dt/taun))
}