NEURON {
    SUFFIX na3
    USEION na READ ena WRITE ina
    RANGE  gbar, ar, sh, shx, mtaufac, htaufac, qa, ina3, thegna, m, h, s, m3h,qinf, thinf
    GLOBAL minf, hinf, mtau, htau, sinf, taus
}

PARAMETER {
    sh   = 0    (mV)        
    shx  = 0    (mV)        
    gbar = 0.010    (mho/cm2)
    mtaufac = 1     (1)     
    htaufac = 1      (1)    
                                
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
    qq   = 10        (mV)
    tq   = -55      (mV)

    thinf  = -62    (mV)        
    qinf  =  6.9    (mV)        

        vhalfs=-60  (mV)        
        a0s=0.0003  (ms)        
        zetas=12    (1)
        gms=0.2     (1)
        smax=10     (ms)
        vvh=-58     (mV) 
        vvs=2       (mV)
        ar=1        (1)     
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
    ina3        (mA/cm2)
    thegna      (mho/cm2)
    minf        hinf 
    m3h     
    mtau (ms)   htau (ms)   
    sinf (ms)   taus (ms)
}
 

STATE { m h s}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h*s
        m3h = m*m*m*h*s
    ina = thegna * (v - ena)
    ina3 = thegna * (v - ena)
    ina=ina3
} 

INITIAL {
    trates(v,ar,sh, shx, mtaufac, htaufac, qa, qinf)
    m=minf  
    h=hinf
    s=sinf
}


FUNCTION alpv(v(mV)) {
         alpv = 1/(1+exp((v-vvh-sh)/vvs))
}
        
FUNCTION alps(v(mV)) {  
  alps = exp(1.e-3*zetas*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bets(v(mV)) {
  bets = exp(1.e-3*zetas*gms*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

LOCAL mexp, hexp, sexp

DERIVATIVE states {   
        trates(v,ar,sh, shx, mtaufac, htaufac, qa, qinf)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        s' = (sinf - s)/taus
}



PROCEDURE trates(vm,a2,sh2, sh3, taufac, taufac2, qa, qinf) {  
        LOCAL  a, b, c, qt
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
    c=alpv(vm)
        sinf = c+a2*(1-c)
        taus = bets(vm)/(a0s*(1+alps(vm)))
        if (taus<smax) {taus=smax}
}




FUNCTION efun(z) {
    if (fabs(z) < 1e-6) {
        efun = 1 - z/2
    }else{
        efun = z/(exp(z) - 1)
    }
}