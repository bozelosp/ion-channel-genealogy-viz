NEURON {
    SUFFIX dIN_kSlow
    USEION k READ ek WRITE ik
    RANGE gmax, ik, n
    RANGE alpha_A, alpha_B, alpha_C, alpha_D, alpha_E
    RANGE beta_A, beta_B, beta_C, beta_D, beta_E
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gmax = 0  (S/cm2)

    alpha_A = 4.61973318e-01
    alpha_B = 8.20458521e-03
    alpha_C = 4.59367292e+00
    alpha_D = -4.20812882e+00
    alpha_E = -1.19678988e+01
    
    beta_A = 9.24268986e-02
    beta_B = -1.35395653e-03
    beta_C = 1.61527685e+00
    beta_D = 2.10656266e+05
    beta_E = 3.32762273e+05
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    ks_ninf
    ks_ntau (ms)
}

STATE {
    n
}

INITIAL {
    ek = -81.5
    rates()
    n = ks_ninf
    
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gmax * pow(n,2) * (v - ek)
}

DERIVATIVE states {
    rates()
    n' =  ( ks_ninf  -  n ) /  ks_ntau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_n, beta_n
    
    

    alpha_n = (alpha_A+alpha_B*v)/(alpha_C+exp((alpha_D+v)/alpha_E))
    beta_n = (beta_A+beta_B*v)/(beta_C+exp((beta_D+v)/beta_E))

    ks_ntau = 1/( alpha_n + beta_n)
    ks_ninf = alpha_n * ks_ntau
    
}

FUNCTION alphabeta(A,B,C,D,E){
    alphabeta = (A + B*v)/(C + exp((v+D)/E))
}
UNITSON