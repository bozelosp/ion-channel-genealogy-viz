UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Nasoma
        USEION na READ ena WRITE ina
        RANGE gna, ina
        GLOBAL minf, hinf, hexp, mtau, htau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 24 (degC)
        dt (ms)
        gna = .0107 (mho/cm2)
        ena = 90 (mV)
}
 
STATE {
        m h 
}
 
ASSIGNED {
        ina (mA/cm2)
        minf 
	mexp 
	hinf 
	hexp
	mtau (ms)
	htau (ms)
}
 
INITIAL {
	rate(v)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE state METHOD cnexp
	ina = gna*m*m*m*h*(v - ena)    
}

DERIVATIVE state {
	rate(v)
	m'=(minf-m)/mtau
	h'=(hinf-h)/htau
}

UNITSOFF
PROCEDURE rate(v(mV)) {  
		      
                      
        LOCAL q10, tinc, alpha, beta
        TABLE minf, hinf, hexp, mtau, htau DEPEND celsius FROM -200 TO 100 WITH 300
		q10 = 3^((celsius - 24)/10)
		tinc = -dt*q10
		alpha = 0.1*vtrap(-(v+38),10)
		beta = 4*exp(-(v+63)/18)
		mtau = 1/(alpha + beta)
		minf = alpha*mtau
		alpha = 0.07*exp(-(v+63)/20)
		beta = 1/(1+exp(-(v+33)/10))
		htau = 1/(alpha + beta)
		hinf = alpha*htau
		hexp = 1-exp(tinc/htau)
}
FUNCTION vtrap(x,y) {	
		if (fabs(x/y) < 1e-6) {
			vtrap = y*(1 - x/y/2)
		}else{
			vtrap = x/(exp(x/y) - 1)
		}
}
UNITSON