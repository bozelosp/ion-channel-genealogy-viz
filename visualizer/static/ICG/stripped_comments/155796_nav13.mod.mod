NEURON {
	SUFFIX nav13
	USEION na READ ena WRITE ina
	RANGE  gbar, ar2, sh, thegna, jina13
	GLOBAL minf, hinf, mtau, htau, sinf, taus,qinf, thinf
	THREADSAFE minf, hinf, mtau, htau, sinf, taus,qinf, thinf
}

PARAMETER {
	sh   = 8	(mV)
	gbar = 1e-1   	(mho/cm2)	
								
	tha  =  -37.5	(mV)		
	qa   = 4.5	(mV)		
	Ra   = 0.4	(/ms)		
	Rb   = 0.135 	(/ms)		

	thi1  = -30	(mV)		
	thi2  = -30 	(mV)		
	qd   = 1.5	(mV)	        
	qg   = 1.5      (mV)
	mmin=0.02	
	hmin=0.5			
	q10=2
	Rg   = 0.01 	(/ms)		
	Rd   = 0.03 	(/ms)		
	qq   = 10        (mV)
	tq   = -50      (mV)

	thinf  = -55 	(mV)		
	qinf  = 4 	(mV)		

        vhalfs=-60	(mV)		
        a0s=0.0003	(ms)		
        zetas=12	(1)
        gms=0.2		(1)
        smax=10		(ms)
        vvh=-58		(mV) 
        vvs=2		(mV)
        ar2=1		(1)		
	ena		(mV)	
	Ena = 55	(mV)            
	celsius
	v 		(mV)
	jina13 		(mA/cm2)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 		
	hinf 		
	mtau 		(ms)	
	htau 		(ms) 	
	sinf 		(ms)	
	taus 		(ms)
}
 

STATE { m h s}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h*s
	ina = thegna * (v - ena)
	jina13 = ina
} 

INITIAL {
	trates(v,ar2,sh)
	m=minf  
	h=hinf
	s=sinf
}


FUNCTION alpv(v(mV)) {
        alpv = 1/(1+exp((v-vvh-sh)/vvs))
}
        
FUNCTION alps(v(mV)) {  
  	alps = exp(1.e-3*zetas*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bets(v(mV)) {
  	bets = exp(1.e-3*zetas*gms*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

LOCAL mexp, hexp, sexp

DERIVATIVE states {   
        trates(v,ar2,sh)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        s' = (sinf - s)/taus
}

PROCEDURE trates(vm,a2,sh2) {  
        LOCAL  a, b, c, qt
        qt=q10^((celsius-24)/10)
	a = trap0(vm,tha+sh2,Ra,qa)
	b = trap0(-vm,-tha-sh2,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = trap0(vm,thi1,Rd,qd) 
	b = trap0(-vm,-thi2,Rg,qg) 
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf)/qinf))
	c=alpv(vm)
        sinf = c+a2*(1-c)
        taus = bets(vm)/(a0s*(1+alps(vm)))
        if (taus<smax) {taus=smax}
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}