INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX sna
    USEION na READ ena WRITE ina     
    GLOBAL minf, hinf, mtau, htau,am,bm,ah,bh     
    RANGE N, reff, eta, gamma, deterministic, gna
    GLOBAL P_am, P_bm, P_ah, P_bh, hinf_curve
    GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf, vshift
    GLOBAL Ra, Rb, Rd, Rg
    GLOBAL vmin, vmax, q10, orig_temp, wflag, tadj
    GLOBAL DONT_VECTORIZE   
}
  
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

PARAMETER {
    v       (mV)
    dt      (ms)
    area

    vshift = -10  (mV)        
                                
    tha  = -35  (mV)        
    qa   = 9    (mV)            
    Ra   = 0.182    (/ms)   
    Rb   = 0.124    (/ms)   

    thi1  = -50 (mV)        
    thi2  = -75 (mV)        
    qi   = 5    (mV)            
    hinf_curve = 1          
                        
                        
    thinf  = -65    (mV)        
    qinf  = 6.2 (mV)        
    Rg   = 0.0091   (/ms)   
    Rd   = 0.024    (/ms)   

    gamma  =  20      (pS)    
    eta     = 50      (1/um2) 
    
    celsius (degC)
    orig_temp = 23 (degC)   
    q10 = 2.3               
    
    deterministic = 0   
    
    vmin = -120 (mV)        
    vmax = 100  (mV)
    DONT_VECTORIZE          
}

ASSIGNED {
    ina         (mA/cm2)
    gna     (pS/um2)
    ena     (mV)

     am     (/ms)
    bm      (/ms)
     ah     (/ms)
    bh      (/ms)
    minf        
     hinf
    mtau (ms)    
     htau (ms)
     tadj

     N 
     reff   (pS/um2)
    scale_dens (pS/um2)
    P_am        
    P_bm
    P_ah
    P_bh
    
    wflag
}
 
STATE { m h                             
        m0h0 m0h1 m1h0 m1h1 m2h0 m2h1 m3h0 m3h1 
m0h0_m1h0  m1h0_m2h0  m2h0_m3h0  m0h1_m1h1  m1h1_m2h1  m2h1_m3h1  
m3h0_m2h0  m2h0_m1h0  m1h0_m0h0  m3h1_m2h1  m2h1_m1h1  m1h1_m0h1  
m0h0_m0h1 m0h1_m0h0 m1h0_m1h1 m1h1_m1h0 m2h0_m2h1 m2h1_m2h0 m3h0_m3h1 m3h1_m3h0 
}


INITIAL { 
    trates(v+vshift) 
    wflag = 1   
    m = minf
    h = hinf
    scale_dens = gamma/area
    reff = gamma*eta
    N   = floor(eta*area + 0.5)   

    m1h0 = floor(3*m*(1-m)*(1-m)*(1-h)*N + 0.5)
    m2h0 = floor(3*m*m*(1-m)*(1-h)*N + 0.5)
    m3h0 = floor(m*m*m*(1-h)*N + 0.5)

    m0h1 = floor((1-m)*(1-m)*(1-m)*h*N + 0.5)
    m1h1 = floor(3*m*(1-m)*(1-m)*h*N + 0.5)
    m2h1 = floor(3*m*m*(1-m)*h*N + 0.5)
    m3h1 = floor(m*m*m*h*N + 0.5)
    
    
    m0h0 = N - (m1h0 + m2h0 + m3h0 + m0h1 + m1h1 + m2h1 + m3h1)

    m0h0_m1h0=0 
    m1h0_m2h0=0 
    m2h0_m3h0=0 
    m0h1_m1h1=0
    m1h1_m2h1=0
    m2h1_m3h1=0 
    m3h0_m2h0=0 
    m2h0_m1h0=0
    m1h0_m0h0=0 
    m3h1_m2h1=0 
    m2h1_m1h1=0 
    m1h1_m0h1=0

    m0h0_m0h1=0 
    m0h1_m0h0=0 
    m1h0_m1h1=0 
    m1h1_m1h0=0 
    m2h0_m2h1=0 
    m2h1_m2h0=0 
    m3h0_m3h1=0 
    m3h1_m3h0=0
}



BREAKPOINT {
    SOLVE states
    if (deterministic) { 
        if (deterministic-1){        
    gna = m*m*m*h*reff*tadj    
    } else {                            
    gna = floor(m*m*m*h* N + 0.5) * scale_dens *tadj}
    } else{                                        
    gna = strap(m3h1) * scale_dens * tadj
    }
    ina = (1e-4) * gna * (v - ena)
} 




PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM
        
    trates(v+vshift)

    
    
    m = m + (1 - exp(-dt/mtau)) * (minf-m)
    h = h + (1 - exp(-dt/htau)) * (hinf-h)
    
    P_am = strap(am*dt)
    P_bm  = strap(bm*dt)
    
    
    ChkProb( 3.0 * P_am)
    ChkProb( 3.0 * P_bm)
    ChkProb( P_bm/(1.0-2.0*P_am) )
    ChkProb( 2.0 * P_bm/(1.0-P_am) )
    
    















    m0h0_m1h0 = BnlDev_RNG(3.0*P_am,m0h0)
    m1h0_m2h0 = BnlDev_RNG(2.0*P_am,m1h0)
    m1h0_m0h0 = BnlDev_RNG(P_bm/(1.0-2.0*P_am), m1h0 - m1h0_m2h0)  
    m2h0_m3h0 = BnlDev_RNG(P_am,m2h0)
    m2h0_m1h0 = BnlDev_RNG(2.0*P_bm/(1.0-P_am), m2h0 - m2h0_m3h0)
    m3h0_m2h0 = BnlDev_RNG(3.0*P_bm,m3h0)
    m0h1_m1h1 = BnlDev_RNG(3.0*P_am, m0h1)
    m1h1_m2h1 = BnlDev_RNG(2.0*P_am, m1h1)
    m1h1_m0h1 = BnlDev_RNG(P_bm/(1.0-2.0*P_am), m1h1 - m1h1_m2h1)
    m2h1_m3h1 = BnlDev_RNG(P_am,m2h1)
    m2h1_m1h1 = BnlDev_RNG(2.0*P_bm/(1.0-P_am), m2h1 - m2h1_m3h1)
    m3h1_m2h1  = BnlDev_RNG(3.0*P_bm,m3h1)

    
    m0h0 = m0h0 - m0h0_m1h0 + m1h0_m0h0
    m1h0 = m1h0 - m1h0_m2h0 - m1h0_m0h0  + m2h0_m1h0 + m0h0_m1h0
    m2h0 = m2h0 - m2h0_m3h0 - m2h0_m1h0  + m3h0_m2h0 + m1h0_m2h0
    m3h0 = m3h0 - m3h0_m2h0 + m2h0_m3h0

    m0h1 = m0h1 - m0h1_m1h1 + m1h1_m0h1
    m1h1 = m1h1 - m1h1_m2h1 - m1h1_m0h1 + m2h1_m1h1 + m0h1_m1h1
    m2h1 = m2h1 - m2h1_m3h1 - m2h1_m1h1 + m3h1_m2h1 + m1h1_m2h1
    m3h1 = m3h1 - m3h1_m2h1 + m2h1_m3h1

    
    P_ah = strap(ah*dt)
    P_bh = strap(bh*dt)
    
    ChkProb(P_ah)
    ChkProb(P_bh)
    
    











    m0h0_m0h1 = BnlDev_RNG(P_ah,m0h0)
    m0h1_m0h0 = BnlDev_RNG(P_bh,m0h1)
    m1h0_m1h1 = BnlDev_RNG(P_ah,m1h0)
    m1h1_m1h0 = BnlDev_RNG(P_bh,m1h1)
    m2h0_m2h1 = BnlDev_RNG(P_ah,m2h0)
    m2h1_m2h0 = BnlDev_RNG(P_bh,m2h1)
    m3h0_m3h1 = BnlDev_RNG(P_ah,m3h0)
    m3h1_m3h0 = BnlDev_RNG(P_bh,m3h1)

    m0h0 = m0h0 - m0h0_m0h1  + m0h1_m0h0
    m1h0 = m1h0 - m1h0_m1h1  + m1h1_m1h0
    m2h0 = m2h0 - m2h0_m2h1  + m2h1_m2h0
    m3h0 = m3h0 - m3h0_m3h1  + m3h1_m3h0

    m0h1 = m0h1 - m0h1_m0h0  + m0h0_m0h1
    m1h1 = m1h1 - m1h1_m1h0  + m1h0_m1h1
    m2h1 = m2h1 - m2h1_m2h0  + m2h0_m2h1
    m3h1 = m3h1 - m3h1_m3h0  + m3h0_m3h1
}




PROCEDURE trates(vm) {     TABLE minf, mtau, hinf, htau, am, bm, ah, bh, tadj
    DEPEND dt,Ra,Rb,Rd,Rg,tha,thi1,thi2,qa,qi,qinf,q10,orig_temp,celsius, hinf_curve
    FROM vmin TO vmax WITH 199
        tadj = q10^((celsius-orig_temp)/10)
    
    
    am = SigmoidRate(vm,tha,Ra,qa)
    am = am * tadj
    bm = SigmoidRate(-vm,-tha,Rb,qa)
    bm = bm * tadj
    mtau = 1/(am+bm)
    minf = am*mtau
    
    
    ah = SigmoidRate(vm,thi1,Rd,qi)
    ah = ah * tadj
    bh = SigmoidRate(-vm,-thi2,Rg,qi)
    bh = bh * tadj
    htau = 1/(ah+bh)
    
    if (hinf_curve == 0) {
        hinf = ah*htau
    }
    else {
        hinf = 1/(1+exp((vm-thinf)/qinf))
        
        ah = hinf/htau
        bh = 1/htau - ah
    }
    
}





FUNCTION SigmoidRate(v,th,a,q) {
    if (fabs(v-th) > 1e-6) {
        SigmoidRate = a * (v - th) / (1 - exp(-(v - th)/q))
    } else {
        SigmoidRate = a * q
    }
}   




FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
    if (wflag){            
VERBATIM
        fprintf (stderr,"sna.mod:strap: negative value for state");
ENDVERBATIM
    wflag = 0}
    } else { 
        strap = x
    }
}



PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
    if (wflag){
VERBATIM
    fprintf(stderr, "sna.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
    wflag =0}
  } }