INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

NEURON {
    SUFFIX skm
    USEION k READ ek WRITE ik
    RANGE gk, gamma, deterministic,reff
    RANGE N,eta
    GLOBAL ninf, ntau,a,b, P_a,P_b
    GLOBAL Rb, Ra
    GLOBAL vmin, vmax, q10, temp, tadj
    GLOBAL DONT_VECTORIZE   
}
   
PARAMETER {
    v           (mV)
    dt      (ms)
    area
    
    gamma = 40       (pS)    
    eta = 0.25       (1/um2)                                 
    tha  = -30  (mV)        
    qa   = 9    (mV)        
    
    Ra   = 0.001    (/ms)       
    Rb   = 0.001    (/ms)       
    
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
    trates(v)
    n = ninf
    scale_dens = gamma/area
    reff = eta*gamma
    N = floor(eta*area + 0.5)
        
    N1 = floor(n * N+ 0.5)
    N0 = N-N1       
    
    n0_n1 = 0
    n1_n0 = 0
}



BREAKPOINT {
  SOLVE states
  
    if (deterministic) { 
        if (deterministic-1){
    gk = n*reff*tadj
    } else {
    gk = floor(n*N + 0.5) * scale_dens*tadj}
    } else{
    gk =  strap(N1) * scale_dens * tadj
    }
    ik = (1e-4) * gk * (v - ek)
} 
       



PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM

    trates(v)
    
    
    
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

PROCEDURE trates(v) {  
                      

        TABLE ninf, ntau, a, b, tadj
    DEPEND  celsius, temp, Ra, Rb, tha, qa,q10
    FROM vmin TO vmax WITH 199

        tadj = q10^((celsius - temp)/10) 
        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
    a = a * tadj                    
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
    b = b * tadj                     
        ntau = 1/(a+b)
        ninf = a/(a+b)
}             



FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skm.mod:strap: negative state");
ENDVERBATIM
    } else {
        strap = x
    }
}



PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
VERBATIM
    fprintf(stderr, "skm.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}