INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nat 
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar, ina
	GLOBAL v1_2m, v1_2h, km, kh
	RANGE minf, hinf, mtau, htau
	GLOBAL vmin, vmax
}

PARAMETER {
	gbar = 0   	(S/cm2)	
	
	v1_2m  = -29.13	(mV)		
	km   = -8.92	(mV)			


	v1_2h  = -47	(mV)		
	kh   = 5	(mV)			

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
	minf 		hinf
	mtau (ms)	htau (ms)
	
}
 

STATE { m h }

INITIAL { 
	trates(v)
	m = minf
	h = hinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	gna = gbar*m*m*m*h
	ina = gna * (v - ena)
} 

DERIVATIVE states {   
        trates(v)      
        
        m' =  (minf-m)/(mtau)
        h' =  (hinf-h)/(htau)
}

PROCEDURE trates(v) {  
                      
        
    TABLE minf,  hinf, mtau, htau
	DEPEND  v1_2m, v1_2h, km, kh
	
	FROM vmin TO vmax WITH 1600

	rates(v)


}


PROCEDURE rates(vm) {  

	
	minf = 1/(1+exp((vm-v1_2m)/km))
	mtau = 0.13 + 3.43/(1+exp((vm+45.35)/5.98))
	if (mtau<1e-7) {mtau=1e-7}

	hinf = 1/(1+exp((vm-v1_2h)/kh))
	htau = 0.36 + exp((vm+20.65)/-10.47)
	if (htau<1e-7) {htau=1e-7}
}