INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ca
	USEION ca READ eca WRITE ica
	RANGE m, h, gcaH, icaH, gbar
	RANGE minf, hinf, mtau, htau
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0.1   	(pS/um2)	
	vshift = 0	(mV)		

	cao  = 2.0	(mM)	        
	cai		(mM)
						
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
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	PI	= (pi) (1)
} 

ASSIGNED {
	ica 		(mA/cm2)
	icaH 		(mA/cm2)
	gcaH		(pS/um2)
	eca		(mV)
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
        SOLVE states METHOD cnexp
        gcaH = gbar*m*m*h
	icaH = (1e-4) * gcaH * (v - eca)
	ica = icaH
} 

LOCAL mexp, hexp








DERIVATIVE states {
        trates(v+vshift)      
        m' =  (minf-m)/mtau
        h' =  (hinf-h)/htau
}

PROCEDURE trates(v) {  
                      
        
        TABLE minf, hinf, mtau, htau 
	DEPEND  celsius, temp
	
	FROM vmin TO vmax WITH 199

	rates(v)





}


PROCEDURE rates(vm) {  
        LOCAL  a, b

        tadj = q10^((celsius - temp)/10)

	a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1)
	b = 0.94*exp((-75-vm)/17)
	
	mtau = 1/tadj/(a+b)
	minf = a/(a+b)

		

	a = 0.000457*exp((-13-vm)/50)
	b = 0.0065/(exp((-vm-15)/28) + 1)

	htau = 1/tadj/(a+b)
	hinf = a/(a+b)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}