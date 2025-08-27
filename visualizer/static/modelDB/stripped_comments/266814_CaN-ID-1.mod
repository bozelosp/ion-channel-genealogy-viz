NEURON {
	SUFFIX canchan
	USEION ca READ eca WRITE ica
	RANGE  gcanbar, ica, pool, bpool
	
}

PARAMETER {
	gcanbar = 100e-6   	(mho/cm2)	
	pool = 1 (1)
	bpool = 1 (1)							
	eca		(mV)            
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ica 		(mA/cm2)
    ican (mA/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        
  ican = gcanbar*m*m*h * (v - eca)
  ica = ican
} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(vm) {  
        LOCAL  a, b

	a = trap0(vm,20,0.1,10)
	b = 0.4*exp(-(vm+25)/18)
	minf = a/(a+b)
	mtau = 1/(a+b)

	a = 0.01*exp(-(vm+50)/10)
	b = 0.1/(1+exp(-(vm+17)/17))
	hinf = a/(a+b)
	htau = 1/(a+b)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}