INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na3
	USEION na READ ena WRITE ina
	RANGE  gbar, ar2
	GLOBAL minf, hinf, mtau, htau, sinf, taus,qinf, thinf
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

        vhalfs=-60	(mV)		
        a0s=0.0003	(ms)		
        zetas=12	(1)
        gms=0.2		(1)
        smax=10		(ms)
        vvh=-58		(mV) 
        vvs=2		(mV)
        ar2=1		(1)		
	ena		(mV)            

	v 		(mV)
	dt		(ms)
	celsius=24	(degC)
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

INITIAL {
trates(v,ar2)
m=minf  
h=hinf
s=sinf
	ina = thegna * (v - ena)
}

BREAKPOINT {
        SOLVE states
        thegna = gbar*m*m*m*h*s
	ina = thegna * (v - ena)
} 

FUNCTION alpv(v(mV)) {
         alpv = 1/(1+exp((v-vvh)/vvs))
}
        
FUNCTION alps(v(mV)) {  
  alps = exp(1.e-3*zetas*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bets(v(mV)) {
  bets = exp(1.e-3*zetas*gms*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
}

LOCAL mexp, hexp, sexp

PROCEDURE states() {   
        trates(v,ar2)      
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        s = s + sexp*(sinf - s)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(vm,a2) {  
        LOCAL  a, b, c, qt
        qt=q10^((celsius-24)/10)
	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)
        mexp = 1 - exp(-dt/mtau)

	a = trap0(vm,thi1,Rd,qd)
	b = trap0(-vm,-thi2,Rg,qg)
	htau =  1/(a+b)/qt
        if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf)/qinf))
        hexp = 1 - exp(-dt/htau)
	c=alpv(vm)
        sinf = c+a2*(1-c)
        taus = bets(vm)/(a0s*(1+alps(vm)))
        if (taus<smax) {taus=smax}
        sexp = (1 - exp(-dt/taus))
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}