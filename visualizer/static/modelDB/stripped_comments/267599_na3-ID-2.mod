NEURON {
	SUFFIX na3
	USEION na READ ena WRITE ina
	RANGE  gbar, ar2,ina
	GLOBAL minf, hinf, mtau, htau, sinf, taus,qinf, thinf
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	
								
	tha  =  -25(mV)		
	qa   = 7.2 (mV)		
	Ra   = 0.4(/ms)		
	Rb   = 0.124 	(/ms)		

	thi1  = -45	(mV)	
	thi2  = -45 	(mV)		
	qd   = 1.5	(mV)	        
	qg   = 1.5      (mV)
	mmin=0.02	
	hmin=0.5			
	q10=2
	Rg   = 0.01 	(/ms)		
	Rd   = .03 	(/ms)		
	qq   = 10        (mV)
	tq   = -55      (mV)

	thinf  = -50	(mV)		
	
       qinf  = 1(mV)		

        vhalfs=-60	(mV)		
        a0s=0.0003	(ms)		
        zetas=12	(1)
        gms=0.2		(1)
        smax=10		(ms)
        vvh=-58		(mV) 
        vvs=2		(mV)
        ar2=1		(1)		
	ena		(mV)            
	celsius
	v 		(mV)
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
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
	sinf (ms)	taus (ms)
}
 

STATE { m h s}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h*s
	ina = thegna * (v - ena)
} 

INITIAL {
	trates(v,ar2)
	m=minf  
	h=hinf
	s=sinf
}


FUNCTION alpv(v(mV)) {
    LOCAL Arg
    Arg=(v-vvh)/vvs
    
    if (Arg<-50) {alpv=1}
    else if (Arg>50) {alpv=0}
    else {alpv=1/(1+exp(Arg))}
}
        
FUNCTION alps(v(mV)) {
    LOCAL Arg
    Arg=1.e-3*zetas*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius))
    
    if (Arg<-50) {alps=0}
    else if (Arg>50) {alps=exp(50)}
    else {alps=exp(Arg)}
}

FUNCTION bets(v(mV)) {
    LOCAL Arg
    Arg=1.e-3*zetas*gms*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius))
    
    if (Arg<-50) {bets=0}
    else if (Arg>50) {bets=exp(50)}
    else {bets=exp(Arg)}
}

LOCAL mexp, hexp, sexp

DERIVATIVE states {   
        trates(v,ar2)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        s' = (sinf - s)/taus
}

PROCEDURE trates(vm,a2) {  
        LOCAL  a, b, c, qt, Arg
        qt=q10^((celsius-24)/10)
	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = trap0(vm,thi1,Rd,qd)
	b = trap0(-vm,-thi2,Rg,qg)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	Arg=(vm-thinf)/qinf
    
    if (Arg<-50) {hinf=1}
    else if (Arg>50) {hinf=0}
    else {hinf=1/(1+exp(Arg))}
	c=alpv(vm)
        sinf = c+a2*(1-c)
        taus = bets(vm)/(a0s*(1+alps(vm)))
        if (taus<smax) {taus=smax}
}

FUNCTION trap0(v,th,a,q) {
	LOCAL Arg
	if (fabs(v-th) > 1e-6) {
	    Arg=-(v-th)/q
        
        if (Arg<-50) {trap0=a*(v-th)}
        else if (Arg>50) {trap0=0}
        else {trap0=a*(v-th)/(1-exp(Arg))}
	} else {
	        trap0 = a * q
 	}
}