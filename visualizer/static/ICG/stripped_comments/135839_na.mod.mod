NEURON {
    THREADSAFE
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 1000   	(pS/um2)	
	vshift = -10	(mV)		
								
	tha  = -35	(mV)		
	qa   = 9	(mV)		
	Ra   = 0.182	(/ms)		
	Rb   = 0.124	(/ms)		

	thi1  = -50	(mV)		
	thi2  = -75	(mV)		
	qi   = 5	(mV)	        
	thinf  = -65	(mV)		
	qinf  = 6.2	(mV)		
	Rg   = 0.0091	(/ms)		
	Rd   = 0.024	(/ms)		

	temp = 23	(degC)		
	q10  = 2.3			


	vmin = -120	(mV)
	vmax = 100	(mV)
}

ASSIGNED {
	v 		(mV)
	celsius		(degC)
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 
STATE { m h }

INITIAL {
    tadj = q10^((celsius - temp)/(10 (degC))) 

	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = tadj*gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 



DERIVATIVE states {   
        trates(v+vshift)      
        m' =  (minf-m)/mtau
        h' =  (hinf-h)/htau
}

PROCEDURE trates(v (mV)) {  
    TABLE minf,  hinf, mtau, htau
    DEPEND celsius, temp, Ra, Rb, Rd, Rg, tha, thi1, thi2, qa, qi, qinf
    FROM vmin TO vmax WITH 199

	rates(v)





}




UNITSOFF
PROCEDURE rates(vm (mV)) {  
    LOCAL  a, b


    a = Ra * qa * efun((tha - vm)/qa)


    b = Rb * qa * efun((vm - tha)/qa)

    tadj = q10^((celsius - temp)/10)

	mtau = 1/tadj/(a+b)
	minf = a/(a+b)

    


    a = Rd * qi * efun((thi1 - vm)/qi)


    b = Rg * qi * efun((vm - thi2)/qi)

    htau = 1/tadj/(a+b)
    hinf = 1/(1+exp((vm-thinf)/qinf))
}
UNITSON



FUNCTION efun(z) {
	if (fabs(z) < 1e-6) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}