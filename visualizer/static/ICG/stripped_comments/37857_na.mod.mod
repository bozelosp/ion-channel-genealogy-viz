UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX na
        USEION na READ ena WRITE ina
        RANGE gnabar, gna, ina
        GLOBAL hinf, minf, htau, mtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  
        dt (ms)
        ena (mV)
        gnabar =  0.07958 (mho/cm2) <0,1e9>
}

STATE {
        m h
}

ASSIGNED {
    ina (mA/cm2) 
    gna (mho/cm2)
    minf hinf
    mtau (ms) htau (ms)
    }

LOCAL mexp, hexp

BREAKPOINT {
	SOLVE states
    
    gna = gnabar*(m^3)*h
    ina = gna*(v - ena)

}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
    h = hinf
}

PROCEDURE states() {  
	trates(v)      
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  
                      

	q10 = 3^((celsius - 22)/10) 


    minf = 1 / (1+exp(-(v + 38) / 7))
    hinf = 1 / (1+exp((v + 65) / 6))

    mtau =  (10 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.04
    htau =  (100 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.6
}

PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE minf, mexp, hinf, hexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    
        
        

	tinc = -dt * q10
	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
	}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON