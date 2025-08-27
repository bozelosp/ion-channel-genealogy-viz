INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX skca
    USEION k READ ek WRITE ik
    USEION ca READ cai        
    RANGE gk, gamma, eta, N, deterministic,reff
    GLOBAL P_a,P_b, ninf, ntau,a,b     
    GLOBAL Ra, Rb, caix
    GLOBAL vmin, vmax, q10, temp, tadj
    GLOBAL DONT_VECTORIZE   
}
   
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

PARAMETER {
    v           (mV)
    dt      (ms)
    area
    
    gamma = 180       (pS)    
    eta = 0.0556         (1/um2)     
    cai         (mM)
    caix = 1    
                                    
    Ra   = 0.01 (/ms)       
    Rb   = 0.02 (/ms)       
        celsius     (degC)
    temp = 23   (degC)      
    q10  = 2.3          
    deterministic = 0   
    vmin = -120 (mV)    
    vmax = 100  (mV)
    DONT_VECTORIZE      
} 

ASSIGNED {
    a       (/ms)
    b       (/ms)
    ik      (mA/cm2)
    gk      (pS/um2)
    ek      (mV)
    ninf        
    ntau (ms)   
    tadj              

    N 
    reff    (pS/um2)
    scale_dens (pS/um2) 
    P_a     
    P_b     
}
 

STATE {
    n         
    N0 N1       
    n0_n1 n1_n0 
}


INITIAL { 
    rates(cai)
    n = ninf
    scale_dens = gamma/area
    reff = eta*gamma
    N = floor(eta*area + 0.5)
    
    N1 = floor(n * N + 0.5)
    N0 = N-N1       
    
    n0_n1 = 0
    n1_n0 = 0
}



BREAKPOINT {
  SOLVE states
   if (deterministic) { 
        if (deterministic-1){
    gk =  n *reff * tadj
    } else { 
    gk = floor(n* N + 0.5) * scale_dens *tadj} 
    } else{                                         
    gk =  strap(N1) * scale_dens * tadj
    }
    ik = (1e-4) * gk * (v - ek)    
}




PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM

    rates(cai)
    
    
    
    n = n + (1 - exp(-dt/ntau)) * (ninf-n)

    P_a = strap(a*dt)
    P_b = strap(b*dt)

    
    ChkProb( P_a)
    ChkProb( P_b)
    
    





    n0_n1 = BnlDev_RNG(P_a, N0)
    n1_n0 = BnlDev_RNG(P_b, N1)

    
    N0    = strap(N0 - n0_n1 + n1_n0)
    N1    = N - N0
}

PROCEDURE rates(cai(mM)) {
        
        tadj = q10^((celsius - temp)/10)
                                  
        a = Ra * cai^caix
    a = a * tadj                                     
        b = Rb
    b = b * tadj                                      
        ntau = 1/(a+b)
        ninf = a/(a+b)
}                                   



FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skca.mod:strap: negative state");
ENDVERBATIM
    } else {
        strap = x
    }
}



PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
VERBATIM
    fprintf(stderr, "skca.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}