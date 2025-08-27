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
      

    SUFFIX NaxSH0_ChannelML
    USEION na READ ena WRITE ina VALENCE 1 ? reversal potential of ion is read, outgoing current is written
            
    RANGE gmax, gion
    
    RANGE minf, mtau
    RANGE hinf, htau
}

PARAMETER { 
      

    gmax = 0.12 (S/cm2) ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
    
    ? Reversal potential of na
    ena (mV)
    ? The outward flow of ion
    ina (mA/cm2)
            
    
    gion (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
         

    gion = gmax*((1*m)^3)*((1*h)^1)
    ina = gion*(v - ena)
                

}



INITIAL {
    ena = 50
        
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
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_m, A_alpha_m, k_alpha_m, d_alpha_m, A_beta_m, k_beta_m, d_beta_m, A_tau_m, k_tau_m, d_tau_m, temp_adj_h, A_alpha_h, k_alpha_h, d_alpha_h, A_beta_h, k_beta_h, d_beta_h, A_tau_h, k_tau_h, d_tau_h, A_inf_h, k_inf_h, d_inf_h
        
    TABLE minf, mtau,hinf, htau
 DEPEND celsius
 FROM -100 TO 100 WITH 2000
    
    
    UNITSOFF
    
    ? There is a Q10 factor which will alter the tau of the gates 
                 

    temp_adj_m = 2^((celsius - 24)/10)     

    temp_adj_h = 2^((celsius - 24)/10)
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_m = 2.880000018
    k_alpha_m = 0.1388888888
    d_alpha_m = -30
     
    
    alpha = A_alpha_m * vtrap((v - d_alpha_m), (1/k_alpha_m))
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_m = 0.892800005
    k_beta_m = -0.1388888888
    d_beta_m = -30
     
    
    beta = A_beta_m * vtrap((v - d_beta_m), (1/k_beta_m))
    
         

    ? Found a generic form of the rate equation for tau, using expression
                    
    
    if (1/( (alpha + beta) * temp_adj_m ) < 0.02 ) {
        tau =  (0.02 * temp_adj_m) 
    } else {
        tau =  1/(alpha + beta)  
    }
    mtau = tau/temp_adj_m
    minf = alpha/(alpha + beta)
          
       
    
    ?     *** Finished rate equations for gate
    
        
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_h = 0.045
    k_alpha_h = 0.6666666667
    d_alpha_h = -45
     
    
    alpha = A_alpha_h * vtrap((v - d_alpha_h), (1/k_alpha_h))
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_h = 0.015
    k_beta_h = -0.6666666667
    d_beta_h = -45
     
    
    beta = A_beta_h * vtrap((v - d_beta_h), (1/k_beta_h))
    
         

    ? Found a generic form of the rate equation for tau, using expression
                    
    
    if (1/( (alpha + beta) * temp_adj_h ) < 0.5 ) {
        tau =  (0.5 * temp_adj_h) 
    } else {
        tau =  1/(alpha + beta)  
    }
    htau = tau/temp_adj_h
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_h = 1
    k_inf_h = 0.25
    d_inf_h = -50
     
    
    inf = A_inf_h / (exp((v - d_inf_h) * k_inf_h) + 1)
    
    hinf = inf
          
       
    
    ?     *** Finished rate equations for gate
    
             

}


? Function to assist with parameterised expressions of type linoid/exp_linear

FUNCTION vtrap(VminV0, B) {
    if (fabs(VminV0/B) < 1e-6) {
    vtrap = (1 + VminV0/B/2)
}else{
    vtrap = (VminV0 / B) /(1 - exp((-1 *VminV0)/B))
    }
}

UNITSON