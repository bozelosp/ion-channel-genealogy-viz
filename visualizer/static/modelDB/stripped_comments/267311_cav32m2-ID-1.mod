INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    THREADSAFE
    SUFFIX it_cav32m2
    USEION ca READ cai,cao WRITE ica
    RANGE pcabar, ica
    GLOBAL m_inf, h_inf, tau_m, tau_h
}

UNITS {
    (molar) = (1/liter)
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (mM)    = (millimolar)
    FARADAY = (faraday) (coulomb)
    R       = (k-mole) (joule/degC)
}

PARAMETER {
    v                   (mV)
    cai                 (mM)
    cao                 (mM)
    celsius             (degC)
    pcabar    = 2e-4    (cm/s)    
}

STATE {
    m h
}

ASSIGNED {
    ica      (mA/cm2)
    m_inf
    tau_m    (ms)
    h_inf
    tau_h    (ms)
}

BREAKPOINT {
    SOLVE castate METHOD cnexp
    ica = pcabar * m * m * h * ghk(v, 1e-6, cao) 
}

DERIVATIVE castate {
    evaluate_fct(v)

    m' = (m_inf - m) / tau_m
    h' = (h_inf - h) / tau_h
}

INITIAL {
    evaluate_fct(v)
    m = m_inf
    h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) {




    UNITSOFF
    m_inf = 1.0 / ( 1.0 + exp(-(v+51.5)/5.8) )
    h_inf = 1.0 / ( 1.0 + exp((v+73.7)/9.1) )
    tau_m = 0.743 + 5.938 / (exp(-(v+70.45)/12.85) + exp((v+70.45)/12.85))
    tau_h = 1.0 + 11.65*exp(-(v+60.0)/233.0)
    UNITSON
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {


    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(z)
    eci = ci*efun(-z)
    ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
        efun = 1 - z/2
    }else{
        efun = z/(exp(z) - 1)
    }
}

UNITSON