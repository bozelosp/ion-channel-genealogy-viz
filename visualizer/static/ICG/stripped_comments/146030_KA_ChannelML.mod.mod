?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file








UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {
      

    SUFFIX KA_ChannelML
    USEION k READ ek WRITE ik VALENCE 1 ? reversal potential of ion is read, outgoing current is written
            
    RANGE gmax, gion
    
    RANGE minf, mtau
    RANGE hinf, htau
}

PARAMETER { 
      

    gmax = 0.0020 (S/cm2) ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
    
    ? Reversal potential of k
    ek (mV)
    ? The outward flow of ion
    ik (mA/cm2)
            
    
    gion (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
         

    gion = gmax*((1*m)^1)*((1*h)^1)
    ik = gion*(v - ek)
                

}



INITIAL {
    ek = -90
        
    rates(v)
    m = minf
                
        
    h = hinf
                
        
    
    
}
    
STATE {
    m
    h
    
}

DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
    
}

PROCEDURE rates(v(mV)) {  
    
    ? Note
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_m, A_alpha_m, k_alpha_m, d_alpha_m, A_beta_m, k_beta_m, d_beta_m, A_tau_m, k_tau_m, d_tau_m, A_inf_m, k_inf_m, d_inf_m, temp_adj_h, A_alpha_h, k_alpha_h, d_alpha_h, A_beta_h, k_beta_h, d_beta_h, A_tau_h, k_tau_h, d_tau_h, A_inf_h, k_inf_h, d_inf_h
        
    TABLE minf, mtau,hinf, htau
 DEPEND celsius
 FROM -100 TO 100 WITH 400
    
    
    UNITSOFF
    
    ? There is a Q10 factor which will alter the tau of the gates 
                 

    temp_adj_m = 3^((celsius - 24)/10)     

    temp_adj_h = 3^((celsius - 24)/10)
    
    ? There is a voltage offset of 0. This will shift the dependency of the rate equations 
    v = v - (0)
    
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_m = 1
    k_alpha_m = 0.1
    d_alpha_m = -45
     
    
    alpha = A_alpha_m * exp((v - d_alpha_m) * k_alpha_m)
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_m = 1
    k_beta_m = 0.075
    d_beta_m = -45
     
    
    beta = A_beta_m * exp((v - d_beta_m) * k_beta_m)
    
         

    ? Found a generic form of the rate equation for tau, using expression
                    tau = beta / (0.04 *(1+alpha))
        
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_m = 1
    k_inf_m = -(0.071428571)
    d_inf_m = 17.5
     
    
    inf = A_inf_m / (exp((v - d_inf_m) * k_inf_m) + 1)
    
    minf = inf
          
       
    
    ?     *** Finished rate equations for gate
    
        
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_h = 1
    k_alpha_h = 0.2
    d_alpha_h = -70
     
    
    alpha = A_alpha_h * exp((v - d_alpha_h) * k_alpha_h)
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_h = 1
    k_beta_h = 0.198
    d_beta_h = -70
     
    
    beta = A_beta_h * exp((v - d_beta_h) * k_beta_h)
    
         

    ? Found a generic form of the rate equation for tau, using expression
                    tau = beta / (0.018 *(1+alpha))
        
    htau = tau/temp_adj_h
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_h = 1
    k_inf_h = (0.166666666)
    d_inf_h = -41.7
     
    
    inf = A_inf_h / (exp((v - d_inf_h) * k_inf_h) + 1)
    
    hinf = inf
          
       
    
    ?     *** Finished rate equations for gate
    
             

}


UNITSON