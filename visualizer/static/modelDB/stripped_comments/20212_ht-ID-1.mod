INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX iar
        USEION h READ eh WRITE ih VALENCE 1
        USEION ca READ cai
        RANGE gbar, h_inf, tau_s, m, shift
        GLOBAL k2, cac, k4, Pc, nca, nexp, ginc, taum
}

UNITS {
        (molar) = (1/liter)
        (mM)    = (millimolar)
        (mA)    = (milliamp)
        (mV)    = (millivolt)
        (msM)   = (ms mM)
}


PARAMETER {
        eh      = -40   (mV)
        celsius = 36    (degC)
        gbar   = 2e-5 (mho/cm2)
        cac     = 0.002 (mM)            
        k2      = 0.0004 (1/ms)         
        Pc      = 0.01                  
        k4      = 0.001 (1/ms)          
        nca     = 4                     
        nexp    = 1                     
        ginc    = 2                     
        taum    = 20    (ms)            
        shift   = 0     (mV)            
}


STATE {
        c1      
        o1      
        o2      
        p0      
        p1      
}


ASSIGNED {
        v       (mV)
        cai     (mM)
        ih      (mA/cm2)
        gh      (mho/cm2)
        h_inf
        tau_s   (ms)
        alpha   (1/ms)
        beta    (1/ms)
        k1ca    (1/ms)
        k3p     (1/ms)
        m
        tadj
}


BREAKPOINT {
        SOLVE ihkin METHOD sparse

        m = o1 + ginc * o2

        ih = gbar * m * (v - eh)
}

KINETIC ihkin {









        evaluate_fct(v,cai)

        ~ c1 <-> o1             (alpha,beta)

        ~ p0 <-> p1             (k1ca,k2)

        ~ o1 <-> o2             (k3p,k4)

        CONSERVE p0 + p1 = 1
        CONSERVE c1 + o1 + o2 = 1
}