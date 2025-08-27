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
	gkdrbar=.003 (mho/cm2)
        vhalfn=13   (mV)
        a0n=0.02      (/ms)
        zetan=-3    (1)
        gmn=0.7  (1)
	nmax=2  (1)
	q10=1
}


NEURON {
	SUFFIX kdrDA
	USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar
	GLOBAL ninf,taun
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkdr
        taun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gkdrbar*n
	ik = gkdr*(v-ek)*DA1(t)*DA2(t)

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { 
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
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