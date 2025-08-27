NEURON {
    SUFFIX hpkj_cn
    NONSPECIFIC_CURRENT i
    RANGE ghbar, eh
    GLOBAL ninf, ntau
    
    RANGE gh, gamma_h
    RANGE Nh, one_over_Nh
    RANGE seed
    
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    
    (pS) = (picosiemens)
    
}

PARAMETER {
    v	 	(mV)
    
    ghbar = .0001	(S/cm2)
    
    eh = -30	(mV)
    
    seed = 5061983 (1)
    gamma_h = 10 (pS)
    
}

ASSIGNED {
    i (mA/cm2)
    ninf
    ntau
    
    
    gh (S/cm2)
    Nh (1)
    one_over_Nh (1)
    
    dt (ms)
    area (um2)
    
    tau_hpkj (ms)
    sigma_hpkj (ms2)
    noise_hpkj
    mu_hpkj
    
}

STATE {
    n
    
    z_hpkj
    
}

INITIAL {
    rates(v)
    n = ninf
    
    Nh = ceil(((1e-8)*area)*(ghbar)/((1e-12)*gamma_h))
    one_over_Nh = 1.0 / Nh
    printf("ihpkj>> the number of channels is %.0f.\n", Nh)
    z_hpkj = 0.
    

}

BREAKPOINT {
    SOLVE states
    gh = ghbar * (n + z_hpkj)
    if (gh < 0) {
        gh = 0
    }
    else if (gh > ghbar) {
        gh = ghbar
    }
    i = gh * (v - eh)
}

PROCEDURE states() {
    rates(v)
    n = n + dt * (ninf - n) / ntau
    
    z_hpkj = z_hpkj*mu_hpkj + noise_hpkj
    
}

PROCEDURE rates(v (mV)) {
    ninf = 1/(1+exp((v+90.1)/9.9))
    ntau = 1000 * (.19 + .72*exp(-((v-(-81.5))/11.9)^2))
    
    mu_hpkj = exp(-dt/ntau)
    sigma_hpkj = one_over_Nh * ninf * (1-ninf)
    noise_hpkj = sqrt(sigma_hpkj * (1-mu_hpkj*mu_hpkj)) * normrand(0,1)
   
}