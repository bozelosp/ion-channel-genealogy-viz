INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na1
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 1.0  	(pS/um2)	
	vshift = -3.5 	(mV)		
								
	tha  = -50	(mV)		
	qa   = 9	(mV)		
	Ra   = 0.182	(/mV/ms)	
	Rb   = 0.124	(/mV/ms)	

	thi1  = -42	(mV)		
	thi2  = -75	(mV)		
	qi   = 5	(mV)	        
	thinf  = -65	(mV)		
	qinf  = 6.2	(mV)		
	Rg   = 0.0091	(/mV/ms)	
	Rd   = 0.024	(/mV/ms)	

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
	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states
        gna = tadj*gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 

LOCAL mexp, hexp 

PROCEDURE states() {   
        trates(v+vshift)      
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(v (mV)) {  
                      
        LOCAL tinc
        TABLE minf, mexp, hinf, hexp
	DEPEND dt, celsius, temp, Ra, Rb, Rd, Rg, tha, thi1, thi2, qa, qi, qinf
	
	FROM vmin TO vmax WITH 199

	rates(v)

        tadj = q10^((celsius - temp)/10(degC))
        tinc = -dt * tadj

        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
}


PROCEDURE rates(vm (mV)) {  
        LOCAL  a, b

	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
	mtau = 1/(a+b)
	minf = a*mtau

		

	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/(a+b)
	hinf = 1/(1+exp((vm-thinf)/qinf))
}


FUNCTION trap0(v (mV), th (mV), a (/mV/ms), q (mV)) (/ms) {
	if (fabs(v/th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}