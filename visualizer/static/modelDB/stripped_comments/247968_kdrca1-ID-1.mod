UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {

	v (mV)
        ek (mV)		
	celsius		(degC)
	gbar=.003 (mho/cm2)
        vhalfn = -15
        a0n=0.02      (/ms)
        zetan=-3    (1)
        gmn=0.7  (1)
	nmax=2  (1)
	qt=1
}


NEURON {
	SUFFIX kdr
	USEION k READ ek WRITE ik
        RANGE gkdr, i, gbar
	RANGE ninf,taun
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
	i  (mA/cm2)
        ninf
        gkdr
        taun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gbar*n
	ik = gkdr*(v-ek)
	i = ik

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*(-3)*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*(-3)*(0.7)*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { 
        LOCAL a
        a = alpn(v)
		if (v < -55 ) {              
		ninf = 0
		} else{
		ninf = 1 / ( 1 + exp( ( vhalfn - v ) / 11 ) )
		
        }
		taun = betn(v)/(qt*(0.02)*(1+a))
	if (taun<nmax) {taun=nmax}
}