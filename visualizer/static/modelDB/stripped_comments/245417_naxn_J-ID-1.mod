NEURON {
	SUFFIX naxj
	USEION na READ ena WRITE ina
	RANGE  gbar, sh
	GLOBAL minf, hinf, mtau, htau 
}

PARAMETER {
	sh   = 0	(mV)
	gbar = 0.010   	(mho/cm2)	
								
	tha  = -105.023 (mV)		
	qa   = 17.7094	(mV)		
	Ra   = 93.8285	(/ms)		
	Rb   = 0.168396	(/ms)		
	qb   = 23.2707	(mV)		

		
	thi  = 17.6769	(mV)		
	qd   = 18.706	(mV)	        
	qg   = 13.3097  (mV)
	mmin=0.002	
	hmin=0.05			
	q10=2
	Rg   = 6.62694 	(/ms)		
	Rd   = 0.000354	(/ms)		

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
	a = trap0(vm,-tha-sh2,Ra,qa) 
	b = Rb * exp(-vm/qb)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = Rd * exp(-vm/qd)
	b = Rg/ (exp(-(v+thi)/qg)+1)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = a/(a+b)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}