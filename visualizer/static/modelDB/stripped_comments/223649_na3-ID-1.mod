NEURON {
	SUFFIX na3
	USEION na READ ena WRITE ina
	RANGE  gbar, i 
	GLOBAL minf, hinf, mtau, htau, qinf, thinf
}

PARAMETER {
	
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
	qq   = 10        (mV)
	tq   = -55      (mV)

	thinf  = -50 	(mV)		
	qinf  = 4 	(mV)		

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
	i    		(mA/cm2)
	thegna		(mho/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
	tha1	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
	i = ina
} 

INITIAL {
	trates(v,ar2)
	m=minf  
	h=hinf
}


LOCAL mexp, hexp

DERIVATIVE states {   
        trates(v,ar2)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(vm,a2) {  
        LOCAL  a, b, qt
		qt = 1.6245
		tha1 = tha 
	a = trap0(vm,tha1,Ra,qa)
	b = trap0(-vm,-tha1,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	
	minf  = 1 / ( 1 + exp( ( - v - 38.43 ) / 7.2 ) )

	a = trap0(vm,thi1,Rd,qd)
	b = trap0(-vm,-thi2,Rg,qg)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	
	hinf  = 1 / ( 1 + exp( ( v + 50 ) / 4 ) )
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}