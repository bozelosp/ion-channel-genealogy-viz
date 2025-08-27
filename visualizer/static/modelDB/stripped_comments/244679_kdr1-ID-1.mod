UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Kdr1
	USEION k WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik, alpha, beta ,pdr,vsh1,vsh2,hdr,pp
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .6 (mho/cm2)
        ek	= -85 (mV)
        pp=6
        pdr=1
        hdr=10
        vsh1 (mV)
        vsh2 (mV)

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
	rates(v,vsh1,vsh2,pdr,hdr,pp)
	m = minf
	h = hinf
}

PROCEDURE states() {  
        rates(v,vsh1,vsh2,pdr,hdr,pp)      
        m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}

PROCEDURE rates(v(mV),vsh1 (mV),vsh2 (mV),pdr,hdr,pp) {  
                      
        LOCAL  q10, tinc, tauh, alpha, beta, gamma, zeta, taum, vshi1,vshi2,corr

        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
                vshi1=vsh1+v
                if(vshi1 == 8) {vshi1 = 8.0001}
                     
        alpha = -0.0047*(vshi1-8)/(exp((vshi1-8)/(-12))-0.9999)
	beta = exp((vshi1+127)/(-30))
     corr=1/(1+exp(-0.4*(v+35)))
	minf = (alpha/(alpha+beta))*corr
    gamma = -0.0047*(vshi1+12)/(exp((vshi1+12)/(-12))-0.9999)
    zeta = exp((vshi1+147)/(-30))
    taum = pdr/(gamma + zeta)
    mexp = 1 - exp(pp*tinc/taum)






                
                vshi2=vsh2+v
        hinf = 1.0 / (1+exp((vshi2+25)/4))
        if(v<-25) {
	tauh = 1200/hdr
	}else{
	tauh = 10/hdr
	}
	hexp = 1 - exp(tinc/tauh)
       
}



 
UNITSON