INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar, vshift
	GLOBAL thm1, thm2, qm1, qm2, thi1, thi2, qi, qinf, thinf
	GLOBAL minf, hinf, mtau, htau
	GLOBAL Am1, Am2, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gbar = 90   	(pS/um2)	
	vshift = 0	(mV)		
								
	thm1  = -45.278153	(mV)		
	thm2  = -17.932028	(mV)		
	Am1   = 0.58733599	(/ms)		
	Am2   = 0.52175946	(/ms)		
	qm1   = 7.8924093	(mV)		
	qm2   = 9.9165402	(mV)		

	thi1  = -28.605225	(mV)		
	thi2  = -40.306515	(mV)		
	qi   = 0.10037405	(mV)	        
	thinf = -54.656584	(mV)		
	qinf  = 0.10001356	(mV)		
	Rg   = 0.011631231	(/ms)		
	Rd   = 0.070431049	(/ms)		

	temp = 23	(degC)		
	q10  = 2.3			

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	rates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states  METHOD cnexp
        gna = gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 

DERIVATIVE states {   
        rates(v+vshift)      
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(vm) {  
        LOCAL  a, b

	a = trap0(vm,thm1,Am1,qm1)
	b = trap0(-vm,-thm2,Am2,qm2)
	mtau = 1/(a+b)
	minf = a*mtau

		

	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/(a+b)
	hinf = 1/(1+exp((vm-thinf)/qinf))
}


FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}