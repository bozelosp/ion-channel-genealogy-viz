NEURON {
	SUFFIX na12mut
	USEION na READ ena WRITE ina
	RANGE  gbar, ar2, thegna
	GLOBAL vhalfs,sh,tha,qa,Ra,Rb,thi1,thi2,qd,qg,mmin,hmin,q10,Rg,qq,Rd,tq,thinf,qinf,vhalfs,a0s,zetas,gms,smax,vvh,vvs
}

PARAMETER {
sh   = 8	(mV)
	gbar = 0.010   	(mho/cm2)	
								
	tha  =  -28.76	(mV)		
	qa   = 5.41	(mV)		
	Ra   = 0.3282 (/ms)		
	Rb   = 0.1 	(/ms)		

	thi1  = -37.651	(mV)		
	thi2  = -30 	(mV)		
	qd   = 0.5	(mV)	        
	qg   = 1.5      (mV)
	mmin=0.02	
	hmin=0.01			
	q10=2
	Rg   = 0.000092 	(/ms)		
	Rd   = .02657 	(/ms)		
	qq   = -65        (mV)    
	tq   = -55      (mV)    

	thinf  = -48.4785 	(mV)		
	qinf  = 7.69	(mV)		

        vhalfs=-20	(mV)		
        a0s=0.00011	(ms)		
        zetas=12	(1)
        gms=0.2		(1)
        smax=10		(ms)
        vvh=-10		(mV) 
        vvs=2		(mV)
        ar2=1		(1)		
	ena		(mV)	
	Ena = 55	(mV)            
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
	minf 		
	hinf 		
	mtau (ms)	
	htau (ms) 	
	sinf (ms)	
	taus (ms)
}
 

STATE { m h s}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h*s
	ina = thegna * (v - Ena)
} 

INITIAL {
	trates(v,ar2,sh)
	m=minf  
	h=hinf
	s=sinf
}


FUNCTION alpv(v) {
         alpv = 1/(1+exp((v-vvh-sh)/vvs))
}
        
FUNCTION alps(v) {  
  alps = exp(1.e-3*zetas*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bets(v) {
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
        if (mtau<mmin) {
		mtau=mmin
		}
	minf = a/(a+b)

	a = trap0(vm,thi1,Rd,qd) 
	b = trap0(-vm,-thi2,Rg,qg) 
	htau =  1/(a+b)/qt
        if (htau<hmin) {
		htau=hmin
		}
	hinf = 1/(1+exp((vm-thinf)/qinf))
	c=alpv(vm)
        sinf = c+a2*(1-c)
        taus = bets(vm)/(a0s*(1+alps(vm)))
        if (taus<smax) {
		taus=smax
		}
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}