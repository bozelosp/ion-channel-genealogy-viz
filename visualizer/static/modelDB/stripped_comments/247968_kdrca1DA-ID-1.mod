UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {

	tone_period = 4000    
	DA_period = 500	
	DA_start = 64000		             
	DA_stop = 96000
	DA_ext1 = 196000
	DA_ext2 = 212000	
	
	DA_t1 = 0.95 
	DA_period2 = 100
	DA_start2 = 36000		   			
	DA_t2 = .8           				

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
	SUFFIX kdrDA
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
	ik = gkdr*(v-ek)*DA1(t)*DA2(t)
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


FUNCTION DA1(t) {
	    if (t >= DA_start && t <= DA_stop){ 									
			if ((t/tone_period-floor(t/tone_period)) >= (1-DA_period/tone_period)) {DA1 = DA_t1}
			else if ((t/tone_period-floor(t/tone_period)) == 0) {DA1 = DA_t1}
			else {DA1 = 1}}
		else if (t >= DA_ext1 && t <= DA_ext2){								
			if ((t/tone_period-floor(t/tone_period)) >= (1-DA_period/tone_period)) {DA1 = DA_t1}
			else if ((t/tone_period-floor(t/tone_period)) == 0) {DA1 = DA_t1}
			else {DA1 = 1}}		
		else  {DA1 = 1}
	}
FUNCTION DA2(t) {
	    if (t >= DA_start2 && t <= DA_stop){
			if((t/tone_period-floor(t/tone_period)) >= (1-DA_period2/tone_period)) {DA2 = DA_t2}
			else if ((t/tone_period-floor(t/tone_period)) == 0) {DA2 = DA_t2}
			else  {DA2 = 1}}
		else  {DA2 = 1}
	}