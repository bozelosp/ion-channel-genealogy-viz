UNITS {
        (molar) = (1/liter)
        (mA) = (milliamp)
        (mV) = (millivolt)
        (mM) = (millimolar)
}

NEURON {
        SUFFIX nmda
        USEION nmda READ enmda WRITE inmda VALENCE 2.0
        NONSPECIFIC_CURRENT i
        RANGE gnmdabar, enmda, g, gnmda, g1, eq, erev
        RANGE A_ap,A_bp,C_ap,C_bp
        GLOBAL pinf, pexp
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        dt (ms)
        g = .0000677254 (mho/cm2)
        erev = -65 (mV)
        gnmdabar = 0.006 (mho/cm2)
        enmda = 0 (mV)
        Mg = 1.8
        A_ap = 0.7
        A_bp = 0.1
        C_ap = 17
        C_bp = 17
}

STATE {
        p
}

ASSIGNED {
        inmda (mA/cm2)
        i  (mA/cm2)
        gnmda (mho/cm2)
        g1 (mho/cm2)
        eq (mV)
        pinf pexp
}

BREAKPOINT {
        SOLVE states
        gnmda = gnmdabar*p
        g1 = g + gnmda
        inmda = gnmda*(v - enmda)
        i = g*(v - erev)
        eq = (g*erev+gnmda*enmda)/g1
}

UNITSOFF

INITIAL {
  rates(v)
  p = pinf
}

PROCEDURE states() {  
        rates(v)      
        p = p + pexp*(pinf-p)
}

PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE pinf, pexp DEPEND dt FROM -100 TO 100 WITH 200
        q10 = 1.0
        tinc = -dt * q10
              
        alpha = A_ap*exp(v/C_ap)
        beta = A_bp*Mg*exp(-v/C_bp)
        sum = alpha + beta
        pinf = alpha/sum
        pexp = 1-exp(tinc*sum)
}

UNITSON