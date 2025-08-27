UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
 
 
}

PARAMETER {
    v       (mV)
    dt      (ms)
    area
    celsius     (degC)
    
    gamma =  10 (pS) 
    eta   =  8  (1/um2) 
    deterministic = 0   
        
        vhalfn=11   (mV)
        vhalfl=-56   (mV)
        a0l=0.05      (/ms)
        a0n=0.05    (/ms)
        zetan=-1.5    (1)
        zetal=3    (1)
        gmn=0.55   (1)
        gml=1   (1)
    lmin=2  (mS)
    nmin=0.1  (mS)
    pw=-1    (1)
    tq=-40
    qq=5
    q10=5
    qtl=1
    ek                                    
    vmin = -120 (mV)    
    vmax = 100  (mV)
    DONT_VECTORIZE          
}


NEURON {
    SUFFIX ska
    USEION k READ ek WRITE ik
    RANGE N, eta, gamma, gka, deterministic,reff
    GLOBAL P_an, P_bn, P_al, P_bl
    GLOBAL ninf, linf, ltau, ntau, lmin
    GLOBAL vmin, vmax, q10
    GLOBAL DONT_VECTORIZE   
}

STATE {
    n l  
    N0L0 N1L0 N0L1 N1L1 
    n0l0_n1l0 n0l0_n0l1 
    n1l0_n1l1 n1l0_n0l0
    n0l1_n1l1 n0l1_n0l0
    n1l1_n0l1 n1l1_n1l0
}

ASSIGNED {
    ik (mA/cm2)
        ninf
        linf      
        ltau  (ms)
        ntau   (ms)
        gka  (pS/um2)  
            an      (/ms)
    bn      (/ms)     al      (/ms)
    bl      (/ms)     reff    (pS/um2)

    N 
    scale_dens (pS/um2) 
    P_an     
    P_bn     
    P_al     
    P_bl     

}

INITIAL {
    rates(v)
    n=ninf
    l=linf
    scale_dens = gamma/area
    reff = eta*gamma
    N = floor(eta*area + 0.5)
    
    N1L1 = floor(n* l* N + 0.5)
    N1L0 = floor(n* (1-l)* N + 0.5)
    N0L1 = floor((1-n)* l* N + 0.5)
    N0L0 = N - N1L1 - N1L0 - N0L1  
    
    n0l0_n1l0 = 0
    n0l0_n0l1 = 0
    n1l0_n1l1 = 0
    n1l0_n0l0 = 0
    n0l1_n1l1 = 0
    n0l1_n0l0 = 0
    n1l1_n0l1 = 0
    n1l1_n1l0 = 0
}


BREAKPOINT {
    SOLVE states 
    if (deterministic) { 
        if (deterministic-1){   
      gka = n * l * reff     
      } else {                                    
      gka = floor(n*l * N + 0.5) * scale_dens}
    } else{                                           
      gka = strap(N1L1) * scale_dens
    }
    ik = gka*(v-ek)*(1e-4)
}




PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM
        
    rates(v)

    
    
    n = n + (1 - exp(-dt/ntau)) * (ninf-n)
    l = l + (1 - exp(-dt/ltau)) * (linf-l)
    
    P_an = strap(an*dt)
    P_bn = strap(bn*dt)
    
    
    ChkProb( P_an)
    ChkProb( P_bn)
    
    
    
    n0l0_n1l0 = BnlDev_RNG(P_an, N0L0)    
    n0l1_n1l1 = BnlDev_RNG(P_an, N0L1)
    n1l1_n0l1 = BnlDev_RNG(P_bn, N1L1)
    n1l0_n0l0 = BnlDev_RNG(P_bn, N1L0)

    
    N0L0 = N0L0 - n0l0_n1l0 + n1l0_n0l0
    N1L0 = N1L0 - n1l0_n0l0 + n0l0_n1l0
    
    N0L1 = N0L1 - n0l1_n1l1 + n1l1_n0l1
    N1L1 = N1L1 - n1l1_n0l1 + n0l1_n1l1

    
    P_al = strap(al*dt)
    P_bl  = strap(bl*dt)
    
    ChkProb(P_al)
    ChkProb(P_bl)
    
    

    n0l0_n0l1 = BnlDev_RNG(P_al,N0L0-n0l0_n1l0)
    n1l0_n1l1 = BnlDev_RNG(P_al,N1L0-n1l0_n0l0)
    n0l1_n0l0 = BnlDev_RNG(P_bl,N0L1-n0l1_n1l1)
    n1l1_n1l0 = BnlDev_RNG(P_bl,N1L1-n1l1_n0l1)

    N0L0 = N0L0 - n0l0_n0l1  + n0l1_n0l0
    N1L0 = N1L0 - n1l0_n1l1  + n1l1_n1l0
    
    N0L1 = N0L1 - n0l1_n0l0  + n0l0_n0l1
    N1L1 = N1L1 - n1l1_n1l0  + n1l0_n1l1
 }

FUNCTION alpn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  alpn = exp(1.e-3*zeta*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  betn = exp(1.e-3*zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}    
PROCEDURE rates(v(mV)) { 
        LOCAL a,qt    
        TABLE ntau, ltau, ninf, linf,al,bl,an,bn
        DEPEND q10, celsius, a0n, nmin, lmin,qtl
        FROM vmin TO vmax WITH 199
         
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1 + a)
        ntau = betn(v)/(qt*a0n*(1+a))
    if (ntau<nmin) {ntau=nmin}
        a = alpl(v)
        linf = 1/(1+ a)
    ltau = 0.26*(v+50)/qtl
    if (ltau<lmin/qtl) {ltau=lmin/qtl}   
    al = linf/ltau
    bl = 1/ltau - al
    an = ninf/ntau
    bn = 1/ntau - an
}



FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skaprox.mod:strap: negative state");
ENDVERBATIM
    } else {
        strap = x
    }
}



PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
VERBATIM
    fprintf(stderr, "skaprox.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}