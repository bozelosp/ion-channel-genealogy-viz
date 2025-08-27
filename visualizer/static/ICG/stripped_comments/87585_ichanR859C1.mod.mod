UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (uF) = (microfarad)
    (molar) = (1/liter)
    (nA) = (nanoamp)
    (mM) = (millimolar)
    (um) = (micron)
    (S) = (siemens)
    FARADAY = 96520 (coul)
    R = 8.3134  (joule/degC)

}

 
NEURON { 
    SUFFIX ichanR859C1 
    
    
    
    USEION na READ ena WRITE ina    
    USEION k READ ek WRITE ik
    RANGE gnat, gkf
    RANGE gnatbar, gkfbar
    RANGE gl, el
    RANGE minf, mtau, hinf, htau, sinf, stau, nfinf, nftau, inat, m, h, s
}

 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

 
PARAMETER {
    
    celsius = 6.3 (degC)
    dt (ms) 
    ena  (mV)
    gnatbar = 0.0 (mho/cm2)   
    ek  (mV)
    gkfbar = 1.0 (mho/cm2)
    gl =0.0 (mho/cm2)    
    el (mV)
}


ASSIGNED {
      
    v (mV) 
    gnat (mho/cm2) 
    gkf (mho/cm2)
    ina (mA/cm2)
    ik (mA/cm2)
    il (mA/cm2)
    minf hinf sinf nfinf
    mtau (ms) htau (ms) stau (ms) nftau (ms)
    mexp hexp sexp nfexp
} 


STATE {
    m h s nf
}
 

BREAKPOINT {
    SOLVE states
    gnat = gnatbar*m*m*m*h*s  
    ina = gnat*(v - ena)
    gkf = gkfbar*nf*nf*nf*nf
    ik = gkf*(v-ek)
    il = gl*(v-el)
}

 
UNITSOFF

 
INITIAL {

    trates(v)
    
    m = minf
    h = hinf
    s = sinf

    nf = nfinf
    
    VERBATIM
    return 0;
    ENDVERBATIM
}


PROCEDURE states() {        
                            
    trates(v)           

    m = m + mexp*(minf-m)
    h = h + hexp*(hinf-h)
    s = s + sexp*(sinf-s)
    nf = nf + nfexp*(nfinf-nf)
    
    VERBATIM
    return 0;
    ENDVERBATIM
}
 

LOCAL q10


PROCEDURE rates(v (mV)) {   
                            

    LOCAL  alpha, beta, sum
    q10 = 3^((celsius - 6.3)/10)
    
    
    
    minf = 1/(1+exp(-(v+21.3)*3.5*0.03937))   	
    mtau = 0.15								 

    
    hinf = 1/(1+exp((v+41.9)/6.7))				
    htau = 23.12*exp(-0.5*((v+77.58)/43.92)^2) 
       
    
    sinf = 1/(1+exp((v+46.0)/6.6))				
    stau = 1000*(190.2*exp(-0.5*((v+90.4)/38.9)^2))

    
    alpha = -0.07*vtrap((v+65-47),-6)
    beta = 0.264/exp((v+65-22)/40)
    sum = alpha+beta        
    nftau = 1/sum      
    nfinf = alpha/sum
}
 

PROCEDURE trates(v (mV)) {  
                            
    LOCAL tinc
 
    TABLE minf, mexp, hinf, hexp, sinf, sexp, nfinf, nfexp, mtau, htau, stau, nftau
        DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
    rates(v)    
                
                

    tinc = -dt * q10
    mexp = 1 - exp(tinc/mtau)
    hexp = 1 - exp(tinc/htau)
    sexp = 1 - exp(tinc/stau)
    nfexp = 1 - exp(tinc/nftau)
}

 
FUNCTION vtrap(x,y) {  

    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{  
        vtrap = x/(exp(x/y) - 1)
    }
}
 

UNITSON