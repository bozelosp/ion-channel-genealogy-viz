UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX M99Ka
        USEION k WRITE ik
        RANGE gkbar,gk,ik,tmfac,afac,bfac,thfac,vms,vhs,minf,hinf,mexp,hexp,thmin

}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 35 (degC)
        dt (ms)
        gkbar = 0.048 (mho/cm2)
        ek = -90 (mV)
        tmfac = 4 (1)   
        thfac = 1 (1)   
        thmin = 2 (ms)  
        afac = 1.5 (1)  
        bfac = 0.825 (1)    
    vms = 0 (mV)    
    vhs = 0 (mV)    
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        minf hinf mexp hexp
}
 
BREAKPOINT {
        SOLVE states
        ik = gkbar*m*h*(v - ek)
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
                      
        LOCAL  alpha, beta, sum, tau
        TABLE minf,mexp,hinf,hexp DEPEND dt,vms,vhs,tmfac,afac,bfac,thfac,thmin FROM -100 TO 100 WITH 200
    
        alpha = exp(-0.038*(afac+1/(1+exp((v-vms+40)/5)))*(v-vms-11))
        beta = exp(-0.038*(bfac+1/(1+exp((v-vms+40)/5)))*(v-vms-11))
        sum = 1 + alpha
        minf = 1/sum
        tau = tmfac*beta/sum
        if (tau < 0.1) {tau=0.1}
        mexp = 1 - exp(-dt/tau)
    
        alpha = exp(0.11*(v-vhs+56))
        hinf = 1/(1+alpha)
        tau = thfac*0.26*(v+50)
        if (tau < thmin) {tau=thmin}
        hexp = 1 - exp(-dt/tau)
}
 
UNITSON