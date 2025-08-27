INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar, vshift
	GLOBAL mslp, mcen, ma, mc, mq1, mq2
	GLOBAL hslp, hcen, ha, hc, hq1, hq2
	GLOBAL minf, hinf, mtau, htau, ina
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gbar = 258.272   	(pS/um2)	
	vshift = 0	(mV)		
								
	mslp = 6.7513		(mV)
	mcen = -31.235		(mV)
	ma = 0.23118		(ms)
	mc = 19.098		(mV)
	mq1 = 1		(mV)
	mq2 = 54.927		(mV)

	hslp  = 2.631		(mV)		
	hcen  = -47.953		(mV)		
	ha   = 14.042		(ms)		
	hc = -73.517		(mV)		
	hq1  = 24.053		(mV)		
	hq2   = 19.627		(mV)	        

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
        SOLVE states  METHOD cnexp
        gna = gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 

DERIVATIVE states {   
        trates(v+vshift)      
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}


PROCEDURE trates(v (mV)) {  
                      
        TABLE minf, mtau , hinf, htau
	DEPEND dt, mcen, mslp, ma, mc, mq1, mq2, hcen, hslp, ha, hc, hq1, hq2
	
	FROM vmin TO vmax WITH 199

UNITSOFF
	rates(v)
UNITSON

}

UNITSOFF

PROCEDURE rates(vm) {  
        LOCAL  a, b

	mtau = 	xtau(vm, ma, mc, mq1, mq2)
	minf =  xinf(vm, mcen, mslp)

		

	htau = xtau(vm, ha, hc, hq1, hq2)
	hinf = xinf(-vm, -hcen, hslp)
}


FUNCTION xinf(v, xcen, xslp) {
	xinf = 1/( 1 + exp(-(v-xcen)/xslp) )
}

FUNCTION xtau(v, xa, xc, xq1, xq2) {
	xtau = xa / (exp(-(v-xc)/xq2) + exp((v-xc)/xq1))
}

UNITSON