NEURON {
    SUFFIX nax
    USEION na READ ena WRITE ina
    RANGE  gbar, sh, shx, mtaufac, htaufac, qa, inax, thegna, m, h, m3h, qinf, thinf
    GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
    sh   = -7   (mV)        
    shx  = -7    (mV)        
    gbar = 0.010    (mho/cm2)   
    mtaufac = 1     (1)     
    htaufac = 1     (1)     
                            
    tha  =  -30 (mV)        
    qa   = 7.2  (mV)        
    Ra   = 0.4  (/ms)       
    Rb   = 0.124    (/ms)       

    thi1  = -45 (mV)        
    thi2  = -45     (mV)        
    qd   = 1.5  (mV)            
    qg   = 1.5      (mV)
    mmin=0.02   
    hmin=0.5            
    q10=2
    Rg   = 0.01     (/ms)       
    Rd   = .03  (/ms)       

    thinf  = -50    (mV)        
    qinf  =  4  (mV)        

    ena     (mV)            
    celsius 
    v       (mV)
}


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

ASSIGNED {
    ina         (mA/cm2)
    inax        (mA/cm2)
    thegna      (mho/cm2)
    minf        hinf
    m3h
    mtau (ms)   htau (ms)   
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
        m3h = m*m*m*h
    ina = thegna * (v - ena)
    inax = thegna * (v - ena)
    ina = inax
} 

INITIAL {
    trates(v,sh, shx, mtaufac, htaufac, qa, qinf)
    m=minf  
    h=hinf
}

DERIVATIVE states {   
        trates(v,sh,shx, mtaufac, htaufac, qa, qinf)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}



PROCEDURE trates(vm,sh2, sh3, taufac, taufac2, qa, qinf) {  
        LOCAL  a, b, qt
        qt=q10^((celsius-24)/10)

    a = Ra * qa * efun((tha+sh2 - vm)/qa)

    b = Rb * qa * efun((vm - tha-sh2)/qa)




    mtau = taufac/(a+b)/qt
        if (mtau<mmin) {mtau=taufac*mmin}
    minf = a/(a+b)
    

    a = Rd * qd * efun((thi1+sh3 - vm)/qd)

    b = Rg * qg * efun((vm - thi2-sh3)/qg)
    


    htau =  taufac2/(a+b)/qt
        if (htau<hmin) {htau=taufac2*hmin}
    hinf = 1/(1+exp((vm-thinf-sh3)/qinf))
}
        


FUNCTION efun(z) {
    if (fabs(z) < 1e-6) {
        efun = 1 - z/2
    }else{
        efun = z/(exp(z) - 1)
    }
}