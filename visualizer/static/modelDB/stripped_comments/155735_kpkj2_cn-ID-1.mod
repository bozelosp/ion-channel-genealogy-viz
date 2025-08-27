NEURON {
    SUFFIX kpkj2_cn
    USEION k READ ek WRITE ik
    RANGE gkbar, ik
    GLOBAL ninf, ntau
    
    RANGE gk, gamma_k
    RANGE Nk, one_over_Nk
    RANGE seed
    
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    
    (S) = (siemens)
    (pS) = (picosiemens)
    
}

PARAMETER {
    v		(mV)
    gkbar = .002	(mho/cm2)
    
    nivh = -24	(mV)
    nik = 20.4
    
    ek
    
    seed = 5061983 (1)
    gamma_k = 10 (pS)
    
}

ASSIGNED {
    ik
    ninf
    ntau		(ms)
    

    gk (S/cm2)
    Nk (1)
    one_over_Nk (1)
    
    dt (ms)
    area (um2)
    
    tau1_kpkj2 (ms) tau2_kpkj2 (ms) tau3_kpkj2 (ms) tau4_kpkj2
    sigma1_kpkj2 (ms2) sigma2_kpkj2 (ms2) sigma3_kpkj2 (ms2) sigma4_kpkj2 (ms2)
    noise1_kpkj2 noise2_kpkj2 noise3_kpkj2 noise4_kpkj2
    mu1_kpkj2 mu2_kpkj2 mu3_kpkj2 mu4_kpkj2
    
    
}

STATE {
    n
    
    z1_kpkj2 z2_kpkj2 z3_kpkj2 z4_kpkj2
    
}

INITIAL {
    rates(v)
    n = ninf
    
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))
    one_over_Nk = 1.0 / Nk
    printf("kpkj2>> the number of channels is %.0f.\n", Nk)
    z1_kpkj2 = 0.
    z2_kpkj2 = 0.
    z3_kpkj2 = 0.
    z4_kpkj2 = 0.
    
}

BREAKPOINT {
    SOLVE states
    gk = gkbar * (n*n*n*n + z1_kpkj2+z2_kpkj2+z3_kpkj2+z4_kpkj2)
    if (gk < 0) {
        gk = 0
    }
    else if (gk > gkbar) {
        gk = gkbar
    }
    ik = gk * (v - ek)
}

PROCEDURE states() {
    rates(v)
    n = n + dt * (ninf - n) / ntau
    
    z1_kpkj2 = z1_kpkj2*mu1_kpkj2 + noise1_kpkj2
    z2_kpkj2 = z2_kpkj2*mu2_kpkj2 + noise2_kpkj2
    z3_kpkj2 = z3_kpkj2*mu3_kpkj2 + noise3_kpkj2
    z4_kpkj2 = z4_kpkj2*mu4_kpkj2 + noise4_kpkj2
    
}

PROCEDURE rates(Vm (mV)) {
    LOCAL v,n4,one_minus_n
    v = Vm + 11	
    ninf = 1/(1+exp(-(v-nivh)/nik))
    ntau = 1000 * ntau_func(v)
    
    
    tau1_kpkj2 = ntau
    tau2_kpkj2 = 0.5 * ntau
    tau3_kpkj2 = 0.3333333 * ntau
    tau4_kpkj2 = 0.25 * ntau
    
    mu1_kpkj2 = exp(-dt/tau1_kpkj2)
    mu2_kpkj2 = exp(-dt/tau2_kpkj2)
    mu3_kpkj2 = exp(-dt/tau3_kpkj2)
    mu4_kpkj2 = exp(-dt/tau4_kpkj2)
    
    n4 = ninf*ninf*ninf*ninf
    one_minus_n = 1. - ninf
    sigma1_kpkj2 = one_over_Nk * 4*n4*ninf*ninf*ninf * one_minus_n
    sigma2_kpkj2 = one_over_Nk * 6*n4*ninf*ninf * one_minus_n*one_minus_n
    sigma3_kpkj2 = one_over_Nk * 4*n4*n * one_minus_n*one_minus_n*one_minus_n
    sigma4_kpkj2 = one_over_Nk * n4 * one_minus_n*one_minus_n*one_minus_n*one_minus_n
    
    noise1_kpkj2 = sqrt(sigma1_kpkj2 * (1-mu1_kpkj2*mu1_kpkj2)) * normrand(0,1)
    noise2_kpkj2 = sqrt(sigma2_kpkj2 * (1-mu2_kpkj2*mu2_kpkj2)) * normrand(0,1)
    noise3_kpkj2 = sqrt(sigma3_kpkj2 * (1-mu3_kpkj2*mu3_kpkj2)) * normrand(0,1)
    noise4_kpkj2 = sqrt(sigma4_kpkj2 * (1-mu4_kpkj2*mu4_kpkj2)) * normrand(0,1)

   
}

FUNCTION ntau_func(v (mV)) {
    if (v < -20) {
	ntau_func = .000688 + 1/(exp((v+64.2)/6.5)+exp((v-141.5)/-34.8))
    } else {
	ntau_func = .00016 + .0008*exp(-.0267 * v)
    }
}