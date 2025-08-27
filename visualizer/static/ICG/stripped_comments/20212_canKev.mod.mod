VERBATIM

extern double nrn_ghk(double, double, double, double);

static const char rcsid[]="$Id: hvaccf.mod,v 1.2 2000/09/27 22:45:41 karchie Exp karchie $";

ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX ca
        USEION ca READ cai,cao WRITE ica
        RANGE pcabar, ica, m_inf, h_inf
}

UNITS {
        (mA)    = (milliamp)
        (mV)    = (millivolt)
        (mM)    = (milli/liter)
        FARADAY = 96480 (coul)
        R       = 8.314 (volt-coul/degC)
}

PARAMETER {
        v                       (mV)
        celsius                 (degC)
        dt                      (ms)
        cai             = 5.e-05(mM)
        cao             = 2.5   (mM)
        pcabar = 1.0                 (cm/s)  
        tauM            = 5     (ms)


        vHalfM          = 3   (mV)
        slopeM          = 8.3    (mV)
        tauH            = 0.8   (ms)


        vHalfH          = -39   (mV)  
        slopeH          = 9.2    (mV)  
        tBase           = 23.5  (degC) 

}

STATE {
        m
        h
}

ASSIGNED {
        ica             (mA/cm2)
        m_inf
        h_inf
}

INITIAL {
        LOCAL tadj

        
        
        
        tadj = 3^((celsius-tBase)/10)   
        tauM = tauM / tadj
        tauH = tauH / tadj

        
        rates(v)
        m = m_inf
        h = h_inf
}

BREAKPOINT {
        SOLVE states


VERBATIM
        ica = pcabar * m * m * h * nrn_ghk((v),(cai),(cao),2);
ENDVERBATIM
}

PROCEDURE states() {
        rates(v)
        m = m + (1-exp(-dt/tauM))*(m_inf-m)
        h = h + (1-exp(-dt/tauH))*(h_inf-h)
}


PROCEDURE rates(v(mV)) {
        m_inf = 1/(1+exp(-(v-vHalfM)/slopeM))
        h_inf = 1/(1+exp((v-vHalfH)/slopeH))
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) {
        LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        e = w / (exp(w)-1)
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else
        
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}