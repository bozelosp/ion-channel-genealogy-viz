UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Kdr
	USEION k WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik, alpha, beta
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .6 (mho/cm2)
        ek	= -85 (mV)

}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar *m*m*h
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {  
        rates(v)      
        m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}

PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, tauh, alpha, beta, gamma, zeta, taum

        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
        alpha = -0.0047*(v-8)/(exp((v-8)/(-12))-0.9999)

	beta = exp((v+127)/(-30))
	minf = alpha/(alpha+beta)
        gamma = -0.0047*(v+12)/(exp((v+12)/(-12))-0.9999)
	zeta = exp((v+147)/(-30))
	taum = 1/(gamma + zeta)
	mexp = 1 - exp(tinc/taum)
                
        hinf = 1.0 / (1+exp((v+25)/4))
        if(v<-25) {
	tauh = 1200
	}else{
	tauh = 10
	}
	hexp = 1 - exp(tinc/tauh)
       
}

 
UNITSON