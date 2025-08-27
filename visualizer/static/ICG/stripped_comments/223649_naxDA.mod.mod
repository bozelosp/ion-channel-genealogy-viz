NEURON {
	SUFFIX naxDA
	USEION na READ ena WRITE ina
	RANGE  gbar
	GLOBAL minf, hinf, mtau, htau,thinf, qinf
}

PARAMETER {
	tone_period = 6000   
	DA_period = 2000
	DA_start = 96000		    
	DA_t1 = -0.3            

	DA_period2 = 200
	DA_start2 = 54001		   
	DA_t2 = -0.9           
	
	gbar = 0.010   	(mho/cm2)	
								
	tha  = -30	(mV)		
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
	celsius		(degC)
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
	tha1
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
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
        LOCAL  a, b, qt
        qt=q10^((celsius-24)/10)
		tha1 = tha + DA1(t)+ DA2(t)
	a = trap0(vm,tha1,Ra,qa)
	b = trap0(-vm,-tha1,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = trap0(vm,thi1,Rd,qd)
	b = trap0(-vm,-thi2,Rg,qg)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf)/qinf))
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	

FUNCTION DA1(t) {
	    if (t > DA_start && (t/tone_period-floor(t/tone_period)) >= (1-DA_period/tone_period)) {DA1 = DA_t1}
		else if (t > DA_start && (t/tone_period-floor(t/tone_period)) == 0) {DA1 = DA_t1}
		else  {DA1 = 0}
	}
FUNCTION DA2(t) {
	    if (t > DA_start2 && (t/tone_period-floor(t/tone_period)) >= (1-DA_period2/tone_period)) {DA2 = DA_t2}
		else if (t > DA_start2 && (t/tone_period-floor(t/tone_period)) == 0) {DA2 = DA_t2}
		else  {DA2 = 0}
	}