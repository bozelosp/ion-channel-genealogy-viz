UNITS {
    (mA) =(milliamp)
    (mV) =(millivolt)
    (uF) = (microfarad)
    (molar) = (1/liter)
    (nA) = (nanoamp)
    (mM) = (millimolar)
    (um) = (micron)
    FARADAY = 96520 (coul)
    R = 8.3134  (joule/degC)
}
 
? interface 
NEURON { 
SUFFIX ichan2vip
USEION nat READ enat WRITE inat VALENCE 1
USEION kf READ ekf WRITE ikf  VALENCE 1
USEION ks READ eks WRITE iks  VALENCE 1
NONSPECIFIC_CURRENT il 
RANGE gnat, gkf, gks
RANGE gnatbar, gkfbar, gksbar
RANGE gl, el
RANGE minf, mtau, hinf, htau, nfinf, nftau, inat, ikf, nsinf, nstau, iks
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
    v (mV) 
    celsius = 6.3 (degC)
    dt (ms) 
    enat  (mV)
    gnatbar (mho/cm2)   
    ekf  (mV)
    gkfbar (mho/cm2)
    eks  (mV)
    gksbar (mho/cm2)
    gl (mho/cm2)    
    el (mV)
}
 
STATE {
    m h nf ns
}
 
ASSIGNED {
         
    gnat (mho/cm2) 
    gkf (mho/cm2)
    gks (mho/cm2)

    inat (mA/cm2)
    ikf (mA/cm2)
    iks (mA/cm2)


    il (mA/cm2)

    minf hinf nfinf nsinf
    mtau (ms) htau (ms) nftau (ms) nstau (ms)
    mexp hexp nfexp nsexp
} 

? currents
BREAKPOINT {
    SOLVE states
    gnat = gnatbar*m*m*m*h  
    inat = gnat*(v - enat)
    gkf = gkfbar*nf*nf*nf*nf
    ikf = gkf*(v-ekf)
    gks = gksbar*ns*ns*ns*ns
    iks = gks*(v-eks)

    il = gl*(v-el)
}
 
UNITSOFF
 
INITIAL {
    trates(v)

    m = minf
    h = hinf
    nf = nfinf
    ns = nsinf

    VERBATIM
    return 0;
    ENDVERBATIM
}

? states
PROCEDURE states() {    
    trates(v)   
    m = m + mexp*(minf-m)
    h = h + hexp*(hinf-h)
    nf = nf + nfexp*(nfinf-nf)
    ns = ns + nsexp*(nsinf-ns)
    VERBATIM
    return 0;
    ENDVERBATIM
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  
                      
    LOCAL  alpha, beta, sum
    
    
    q10 = 1     
    
    
    alpha = -0.3*vtrap((v+41.5),-5)
    beta  = 0.3*vtrap((v+25),5)
    sum   = alpha+beta        
    mtau  = 1/sum      
    minf  = alpha/sum
    
    
    alpha = 0.23/exp((v+65)/20)
    beta  = 3.33/(1+exp((v+37.5)/-10))
    sum   = alpha+beta
    htau  = 1/sum 
    hinf  = alpha/sum 
    
    
    alpha = -0.88*vtrap((v+30),-6)
    beta  = 0.1056/exp((v+55)/40)
    sum   = alpha+beta        
    nstau = 1/sum      
    nsinf = alpha/sum
    
    
    alpha = -0.065*vtrap((v+18),-6)
    beta  = 0.264/exp((v+28)/40)
    sum   = alpha+beta        
    nftau = 1/sum      
    nfinf = alpha/sum
}
 
PROCEDURE trates(v) {  
                       
    LOCAL tinc
    TABLE minf, mexp, hinf, hexp, nfinf, nfexp, nsinf, nsexp, mtau, htau, nftau, nstau
    DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
    rates(v)  
              
              

    tinc  = -dt * q10
    mexp  = 1 - exp(tinc/mtau)
    hexp  = 1 - exp(tinc/htau)
    nfexp = 1 - exp(tinc/nftau)
    nsexp = 1 - exp(tinc/nstau)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON