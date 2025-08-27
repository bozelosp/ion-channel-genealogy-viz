NEURON {
	SUFFIX kad
	USEION k READ ek WRITE ik
        RANGE gkabar,gka,ik
        GLOBAL ninf,linf,taul,taun,lmin
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {    

	gkabar = 0      (mho/cm2)  
        vhalfn = -1     (mV)       
        vhalfl = -56    (mV)       
       a0n = 0.1       (/ms)      
       
        zetan = -1.8    (1)        
        zetal = 3       (1)        
        gmn   = 0.39    (1)        
        gml   = 1       (1)
	lmin  = 2       (ms)

	nmin  = 0.2     (ms)	
	pw    = -1      (1)
	tq    = -40     (mV)
	qq    = 5       (mV)
	q10   = 5                
}


ASSIGNED {    
	v               (mV)
        ek              (mV)
	celsius  	(degC)
	ik              (mA/cm2)
        ninf
        linf      
        taul            (ms)
        taun            (ms)
        gka             (mho/cm2)
}


STATE {       
	n l
}


LOCAL qt

INITIAL {    
        qt = q10^((celsius-24)/10(degC))       
	rates(v)
	n=ninf
	l=linf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkabar*n*l*(v+70)

}


DERIVATIVE states {     
        rates(v)          
        n' = (ninf - n)/taun
        l' = (linf - l)/taul
}



PROCEDURE rates(v (mV)) {		 
	LOCAL a

        a = alpn(v)
        ninf = 1/(1 + a)		 
        taun = betn(v)/(qt*a0n*(1+a))	 
	if (taun<nmin) {taun=nmin}	 

        a = alpl(v)
        linf = 1/(1+ a)                  
	taul = 0.26(ms/mV)*(v+50)               
	if (taul<lmin) {taul=lmin}       
}


FUNCTION alpn(v(mV)) { LOCAL zeta
  zeta = zetan+pw/(1+exp((v-tq)/qq))
UNITSOFF
  alpn = exp(1.e-3*zeta*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION betn(v(mV)) { LOCAL zeta
  zeta = zetan+pw/(1+exp((v-tq)/qq))
UNITSOFF
  betn = exp(1.e-3*zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION alpl(v(mV)) {
UNITSOFF
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION betl(v(mV)) {
UNITSOFF
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}