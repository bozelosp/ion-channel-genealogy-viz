UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
        ek (mV)		
	celsius		(degC)
	gkdrbar=.003 (mho/cm2)
       vhalfn=13 (mV)
        a0n=0.02      (/ms)
        zetan=-3    (1)
        gmn=0.7  (1)

	nmax=0.10  (ms) 
        tau (ms)
	q10=1
}


NEURON {
	SUFFIX kdr
	USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar,tau,ik
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        taun (ms)
        gkdr (mho/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
        tau = 0.7*betn(v)/(0.02*(1+ exp(-339.44*(v-13)/(8.315*(273.16+celsius))) ))
	gkdr = gkdrbar*n
	ik = gkdr*(v-ek)

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-13)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { 
        LOCAL a,c,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1+a)
 c = exp(-339.44*(v-13)/(8.315*(273.16+celsius))) 
       taun = 0.7*betn(v)/(qt*a0n*(1+c))

        if(taun<nmax) {taun=nmax}
   
}