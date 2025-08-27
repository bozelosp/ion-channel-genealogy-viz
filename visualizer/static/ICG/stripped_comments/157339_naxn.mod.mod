NEURON {
	SUFFIX nax
	USEION na READ ena WRITE ina
	RANGE  gbar, sh
	GLOBAL minf, hinf, mtau, htau,thinf, qinf, ena
}

PARAMETER {
	sh   = 0	(mV)
	gbar = 0.010   	(mho/cm2)	
								
	tha  =  -30	(mV)		
	qa   = 7.2	(mV)		
	Ra   = 0.4	(/ms)		
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

	thinf  = -50 	(mV)		
	qinf  = 4 	(mV)		

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
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
} 

INITIAL {
	trates(v,sh)
	m=minf  
	h=hinf
}

DERIVATIVE states {   
        trates(v,sh)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(vm,sh2) {  
        LOCAL  a, b, qt
        qt=q10^((celsius-24)/10)
	a = trap0(vm,tha+sh2,Ra,qa)
	b = trap0(-vm,-tha-sh2,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = trap0(vm,thi1+sh2,Rd,qd)
	b = trap0(-vm,-thi2-sh2,Rg,qg)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf-sh2)/qinf))
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}