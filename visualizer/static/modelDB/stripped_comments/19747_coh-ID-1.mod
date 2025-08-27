DEFINE SSIZE 501
DEFINE RSIZE 500

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS COH
    RANGE T, ntot, Ttot, spike
    GLOBAL KC, KO
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mM) = (milli/liter)
}

PARAMETER {
    dt    (ms)
    rseed = 0   
    pv0 = 0.5  (1)  
    n0 = 1    (1)   
    ke = 8 (/mM /ms)   
    kd = 0.0002  (/ms)  
    Camp = 0.1 (mM) 
    Cres = 0 (mM)   
    Cdur = 1 (ms)   
    Cnamp = 0.01 (mM)   
    Cnres = 0 (mM)  
    Cndur = 2 (ms)  
    Tamp = 1 (mM)   
    Tdur = 1  (ms)  
}

ASSIGNED {
    spike[SSIZE]    (ms)    
    index           
    tspike      (ms)    
    trel[RSIZE] (ms)    
    ntot        (1) 
    Ttot        (mM)    
    km      (/ms)   
    KO[2]       (/mM /ms) 
    KC[2]       (/ms)     
    inf[2] tau[2] fac[2]
}

STATE {
    n[RSIZE]    (1) 
    R   (1) 
    RO[2]   (1) 
    T[RSIZE]    (mM)    
    C       (mM)    
    Cn      (mM)    
    pv[RSIZE]   (1) 
    Ccount (1)  
    Cncount (1)  
    Tcnt[RSIZE] (1)  
}

INITIAL {
    index = 0
    FROM i = 0 TO RSIZE-1 {
      n[i] = n0
      T[i] = 0
      trel[i] = 0
      pv[i] = pv0
    }
    R = 0
    C = Cres
    Cn = Cnres
    ntot = 0
    Ttot = 0
    set_seed(rseed)
    km = n0 * kd        
    KO[0] = 150 
    KO[1] = 1 
    KC[0] = 30 
    KC[1] = 0.1 
}

BREAKPOINT {
    SOLVE release
}

PROCEDURE release() {

  if (index < SSIZE && t>=spike[index] && C < Camp) {    
    C = Camp    
    index=index+1 
    tspike = t    
    Ccount = Cdur / dt  
    Ttot = 0
  }

  prel()    

  ntot = 0
  FROM i = 0 TO RSIZE-1 {

    if (unirand() < dt*km) {n[i] = n[i]+1}  
    if (n[i] > 0 && unirand() < dt*kd) {n[i] = n[i]-1}  
    if (Cn > 0) {
      if (unirand() < dt*ke*Cn) {n[i] = n[i]+1} 
    }
    
    onerel(i)   

    if (T[i] > 0) {
      Tcnt[i] = Tcnt[i] - 1
      if (Tcnt[i] < 0) {T[i] = 0}    
    }

    ntot = ntot + n[i] 
  }

  if (C > Cres) {
    Ccount = Ccount - 1
    if (Ccount < 0) {
      C = Cres  
      Cn = Cnamp  
      Cncount = Cndur / dt  
    }
  }
  
  if (Cn > Cnres) {
    Cncount = Cncount - 1 
    if (Cncount < 0) {Cn = Cnres}
  }
    
  VERBATIM
  return 0;
  ENDVERBATIM
}


PROCEDURE onerel(i) {   
  if (trel[i] < tspike) {
    
    if (n[i] > 0 && unirand() < R*pv[i]*n[i]*dt*40) { 
      n[i] = n[i] - 1     
      T[i] = Tamp     
      Tcnt[i] = Tdur / dt     
      trel[i] = t     
      Ttot = Ttot + 1     
    }
  }
}


PROCEDURE prel() {  
  rates(C)
  R = 1
  FROM i=0 TO 1 {
    RO[i] = RO[i] + fac[i]*(inf[i] - RO[i])
    R = R * RO[i]
  }
}

 
PROCEDURE rates(C) {LOCAL a, b  
        TABLE inf, fac, tau DEPEND dt FROM 0 TO 0.2 WITH 200
        FROM j=0 TO 1 {
            a = KO[j] * C
            b = KC[j]
            tau[j] = 1/(a + b)
            inf[j] = a/(a + b)
            fac[j] = (1 - exp(-dt/tau[j]))
        }
}


FUNCTION unirand() {    
        return(scop_random())
}