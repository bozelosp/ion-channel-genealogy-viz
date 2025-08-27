INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS AMPA
    POINTER C
    RANGE C0, C1, C2, D, O1, O2
    RANGE g, gmax, rb
    GLOBAL Erev
    GLOBAL Rb, Ru1, Ru2, Rd, Rr, Ro1, Rc1, Ro2, Rc2
    GLOBAL vmin, vmax
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (umho) = (micromho)
    (mM) = (milli/liter)
    (uM) = (micro/liter)
}

PARAMETER {

    Erev    = 7    (mV) 
    gmax    = 500  (pS)  
    vmin = -120 (mV)
    vmax = 100  (mV)
    


    Rb  = 13   (/mM /ms)
                
    Ru1 = 0.3  (/ms)    
    Ru2 = 200  (/ms)    
    Rd  = 30.0   (/ms)  
    Rr  = 0.02 (/ms)    
    Ro1 = 100    (/ms)  
    Rc1 = 2    (/ms)  
    Ro2 = 2    (/ms)   
    Rc2 = 0.25    (/ms) 
}

ASSIGNED {
    v       (mV)    
    i       (nA)    
    g       (pS)    
    C       (mM)    
    rb      (/ms)   
}

STATE {
    
    C0      
    C1      
    C2      
    D       
    O1      
    O2      
}

INITIAL {
    C0=1
    C1=0
    C2=0
    D=0
    O1=0
    O2=0
}

BREAKPOINT {
    SOLVE kstates METHOD sparse

    g = gmax * (O1 + O2)
    i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
    
    rb = Rb * C 

    ~ C0  <-> C1    (rb,Ru1)
    ~ C1 <-> C2 (rb,Ru2)
    ~ C2 <-> D  (Rd,Rr)
    ~ C2 <-> O1 (Ro1,Rc1)
    ~ C2 <-> O2 (Ro2,Rc2)

    CONSERVE C0+C1+C2+D+O1+O2 = 1
}