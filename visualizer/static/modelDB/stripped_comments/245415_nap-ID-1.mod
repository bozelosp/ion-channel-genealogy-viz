INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE m, gna, gbar, ina
	GLOBAL v1_2m, km
	RANGE minf, mtau
	GLOBAL vmin, vmax
}

PARAMETER {
	gbar = 0   	(S/cm2)	
	
	v1_2m  = -48.77	(mV)		
	km   = -3.68	(mV)			



	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 1000	(mV)
	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(S/cm2)
	ena		(mV)
	minf 		
	mtau 		(ms)	
	
}
 

STATE {m}

INITIAL { 
	trates(v)
	m = minf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	gna = gbar*m
	ina = gna * (v - ena)
} 


DERIVATIVE states {   
        trates(v)      
        
        m' =  (minf-m)/(mtau)
}

PROCEDURE trates(v) {  
                      
        
    TABLE minf, mtau
	DEPEND  v1_2m, km
	
	FROM vmin TO vmax WITH 1600

	rates(v)


}


PROCEDURE rates(vm) {  

	
	minf = 1/(1+exp((vm-v1_2m)/km))
	mtau = 1
}