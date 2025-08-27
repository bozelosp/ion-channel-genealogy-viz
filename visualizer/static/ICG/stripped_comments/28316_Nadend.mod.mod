UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Nadend
        USEION na READ ena WRITE ina
        NONSPECIFIC_CURRENT il
        RANGE gnadend, gl, el, ina
        GLOBAL minf, hinf, hexp, mtau, htau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 24 (degC)
        dt (ms)
        gnadend = .0117 (mho/cm2)
        
        gl = .00005 (mho/cm2)
        el = -70 (mV)
}
 
STATE {
        m h 
}
 
ASSIGNED {
        ena (mV)
        ina (mA/cm2)
        il (mA/cm2)
        minf 
	mexp 
	hinf 
	hexp
	mtau (ms)
	htau (ms)
}
 
INITIAL {
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states
	ina = gnadend*minf*minf*minf*h*(v - ena)    
        il = gl*(v - el)
}

PROCEDURE states() {	
	evaluate_fct(v)
	h = h + hexp*(hinf - h)
	VERBATIM
	return 0;
	ENDVERBATIM 
}
UNITSOFF

PROCEDURE evaluate_fct(v(mV)) {  
		      
                      
        LOCAL q10, tinc, alpha, beta
        TABLE minf, hinf, hexp, mtau, htau DEPEND dt, celsius FROM -200 TO 
100 WITH 300
		q10 = 3^((celsius - 24)/10)
		tinc = -dt*q10
		alpha = 0.1*vtrap(-(v+45),10)
		beta = 4*exp(-(v+70)/18)
		mtau = 1/(alpha + beta)
		minf = alpha*mtau
		alpha = 0.07*Exp(-(v+70)/20)
		beta = 1/(1+Exp(-(v+40)/10))
		htau = 1/(alpha + beta)
		hinf = alpha*htau
		hexp = 1-Exp(tinc/htau)
}
FUNCTION vtrap(x,y) {	
		if (fabs(x/y) < 1e-6) {
			vtrap = y*(1 - x/y/2)
		}else{
			vtrap = x/(Exp(x/y) - 1)
		}
}
FUNCTION Exp(x) {
		if (x < -100) {
			Exp = 0
		}else{
			Exp = exp(x)
		}
}
UNITSON