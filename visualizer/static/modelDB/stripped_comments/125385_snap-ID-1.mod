INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX snap
    USEION na READ ena WRITE ina     
    GLOBAL minf, mtau,am,bm     
    RANGE N, reff, eta, gamma, deterministic, gna
    GLOBAL P_am, P_bm
    GLOBAL tha, qa, vshift
    GLOBAL Ra, Rb
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
    minf 
    mtau (ms)
     tadj

     N 
     reff   (pS/um2)
    scale_dens (pS/um2)
    P_am        
    P_bm
    
    wflag
}
 
STATE { m                             
    m00 m1 m2 m3 
    m00_m1  m1_m2  m2_m3 
    m3_m2  m2_m1  m1_m00 
}



INITIAL { 
    trates(v+vshift) 
    wflag = 1   
    m = minf
    scale_dens = gamma/area
    reff = gamma*eta
    N   = floor(eta*area + 0.5)   

    m1 = floor(3*m*(1-m)*(1-m)*N + 0.5)
    m2 = floor(3*m*m*(1-m)*N + 0.5)
    m3 = floor(m*m*m*N + 0.5)

    
    m00 = N - (m1 + m2 + m3)


    m00_m1 = 0  
    m1_m2 = 0
    m2_m3 = 0
    m3_m2 = 0 
    m2_m1 = 0 
    m1_m00 = 0  
}



BREAKPOINT {
    SOLVE states
    if (deterministic) { 
        if (deterministic-1){        
    gna = m*m*m*reff*tadj    
    } else {                            
    gna = floor(m*m*m*N + 0.5) * scale_dens *tadj}
    } else{                                        
    gna = strap(m3) * scale_dens * tadj
    }
    ina = (1e-4) * gna * (v - ena)
} 




PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM     
    trates(v+vshift)

    
    
    m = m + (1 - exp(-dt/mtau)) * (minf-m)
    
    P_am = strap(am*dt)
    P_bm  = strap(bm*dt)
    
    
    ChkProb( 3.0 * P_am)
    ChkProb( 3.0 * P_bm)
    ChkProb( P_bm/(1.0-2.0*P_am) )
    ChkProb( 2.0 * P_bm/(1.0-P_am) )

    m00_m1 = BnlDev_RNG(3.0*P_am,m00)
    m1_m2 = BnlDev_RNG(2.0*P_am,m1)
    m1_m00 = BnlDev_RNG(P_bm/(1.0-2.0*P_am), m1 - m1_m2)  
    m2_m3 = BnlDev_RNG(P_am,m2)
    m2_m1 = BnlDev_RNG(2.0*P_bm/(1.0-P_am), m2 - m2_m3)
    m3_m2 = BnlDev_RNG(3.0*P_bm,m3)
    
    
    m00 = m00 - m00_m1 + m1_m00
    m1 = m1 - m1_m2 - m1_m00  + m2_m1 + m00_m1
    m2 = m2 - m2_m3 - m2_m1  + m3_m2 + m1_m2
    m3 = m3 - m3_m2 + m2_m3
}




PROCEDURE trates(vm) {     
    TABLE minf, mtau, am, bm, tadj
    DEPEND dt,Ra,Rb,tha,qa,q10,orig_temp,celsius
    FROM vmin TO vmax WITH 199
        tadj = q10^((celsius-orig_temp)/10)
    
    
    am = SigmoidRate(vm,tha,Ra,qa)
    am = am * tadj
    bm = SigmoidRate(-vm,-tha,Rb,qa)
    bm = bm * tadj
    mtau = 1/(am+bm)
    minf = am*mtau
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
        fprintf (stderr,"snap.mod:strap: negative value for state");
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
    fprintf(stderr, "snap.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
    wflag =0}
  } }