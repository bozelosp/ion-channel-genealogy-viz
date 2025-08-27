UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX klt
        USEION k READ ek WRITE ik
        RANGE gkltbar, gklt, ik
        GLOBAL winf, zinf, wtau, ztau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  
        dt (ms)
        ek = -77 (mV)
        gkltbar = 0.01592 (mho/cm2) <0,1e9>
        zss = 0.5   <0,1>   
}

STATE {
        w z
}

ASSIGNED {
    ik (mA/cm2) 
    gklt (mho/cm2)
    winf zinf
    wtau (ms) ztau (ms)
    }

LOCAL wexp, zexp

BREAKPOINT {
	SOLVE states
    
	gklt = gkltbar*(w^4)*z
    ik = gklt*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    w = winf
    z = zinf
}

PROCEDURE states() {  
	trates(v)      
	w = w + wexp*(winf-w)
	z = z + zexp*(zinf-z)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  
                      

	q10 = 3^((celsius - 22)/10) 

    winf = (1 / (1 + exp(-(v + 48) / 6)))^0.25
    zinf = zss + ((1-zss) / (1 + exp((v + 71) / 10)))

    wtau =  (100 / (6*exp((v+60) / 6) + 16*exp(-(v+60) / 45))) + 1.5
    ztau =  (1000 / (exp((v+60) / 20) + exp(-(v+60) / 8))) + 50
}

PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE winf, wexp, zinf, zexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    
        
        

	tinc = -dt * q10
	wexp = 1 - exp(tinc/wtau)
	zexp = 1 - exp(tinc/ztau)
	}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON