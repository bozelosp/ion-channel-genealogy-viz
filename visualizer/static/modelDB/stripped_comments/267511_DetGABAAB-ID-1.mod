NEURON {
    THREADSAFE

    POINT_PROCESS DetGABAAB
    RANGE tau_r_GABAA, tau_d_GABAA, tau_r_GABAB, tau_d_GABAB
    RANGE Use, u, Dep, Fac, u0, GABAB_ratio
    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB
    NONSPECIFIC_CURRENT i
    RANGE synapseID, verboseLevel
}


PARAMETER {
    tau_r_GABAA  = 0.2   (ms) 
    tau_d_GABAA  = 8     (ms) 
    tau_r_GABAB  = 3.5   (ms) 
    tau_d_GABAB  = 260.9 (ms) 
    Use          = 1.0   (1)  
    Dep          = 100   (ms) 
    Fac          = 10    (ms) 
    e_GABAA      = -80   (mV) 
    e_GABAB      = -97   (mV) 
    gmax         = .001  (uS) 
    u0           = 0          
    synapseID    = 0
    verboseLevel = 0
    GABAB_ratio  = 0     (1)  
}


ASSIGNED {
    v (mV)
    i (nA)
    i_GABAA (nA)
    i_GABAB (nA)
    g_GABAA (uS)
    g_GABAB (uS)
    g (uS)
    factor_GABAA
    factor_GABAB
}

STATE {
    A_GABAA       
    B_GABAA       
    A_GABAB       
    B_GABAB       
}


INITIAL{
    LOCAL tp_GABAA, tp_GABAB

    A_GABAA = 0
    B_GABAA = 0

    A_GABAB = 0
    B_GABAB = 0

    tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA) 
    tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) 

    factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) 
    factor_GABAA = 1/factor_GABAA

    factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) 
    factor_GABAB = 1/factor_GABAB
}


BREAKPOINT {
    SOLVE state METHOD cnexp
    g_GABAA = gmax*(B_GABAA-A_GABAA) 
    g_GABAB = gmax*(B_GABAB-A_GABAB) 
    g = g_GABAA + g_GABAB
    i_GABAA = g_GABAA*(v-e_GABAA) 
    i_GABAB = g_GABAB*(v-e_GABAB) 
    i = i_GABAA + i_GABAB
}


DERIVATIVE state{
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
    A_GABAB' = -A_GABAB/tau_r_GABAB
    B_GABAB' = -B_GABAB/tau_d_GABAB
}


NET_RECEIVE (weight,weight_GABAA, weight_GABAB, R, Pr, u, tsyn (ms)){
    LOCAL result
    weight_GABAA = weight
    weight_GABAB = weight * GABAB_ratio

    INITIAL{
            R=1
            u=u0
            tsyn=t
    }

    
    if (Fac > 0) {
        u = u*exp(-(t - tsyn)/Fac) 
    } else {
        u = Use
    }
    if(Fac > 0){
        u = u + Use*(1-u) 
    }

    R  = 1 - (1-R) * exp(-(t-tsyn)/Dep) 
                                        
    Pr  = u * R                         
    R  = R - u * R                      

    if( verboseLevel > 0 ) {
        printf("Synapse %f at time %g
    }

    tsyn = t

    A_GABAA = A_GABAA + Pr*weight_GABAA*factor_GABAA
    B_GABAA = B_GABAA + Pr*weight_GABAA*factor_GABAA
    A_GABAB = A_GABAB + Pr*weight_GABAB*factor_GABAB
    B_GABAB = B_GABAB + Pr*weight_GABAB*factor_GABAB

    if( verboseLevel > 0 ) {
        printf( " vals %g %g %g %g\n", A_GABAA, weight_GABAA, factor_GABAA, weight )
    }
}


FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}