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
      

    SUFFIX Kdr_ChannelML
    USEION k READ ek WRITE ik VALENCE 1 ? reversal potential of ion is read, outgoing current is written
            
    RANGE gmax, gion
    
    RANGE minf, mtau
}

PARAMETER { 
      

    gmax = 0.036 (S/cm2) ? default value, should be overwritten when conductance placed on cell
    
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
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
         

    gion = gmax*((1*m)^1)
    ik = gion*(v - ek)
                

}



INITIAL {
    ek = -77.0
        
    rates(v)
    m = minf
                
        
    
    
}
    
STATE {
    m
    
}

DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    
}

PROCEDURE rates(v(mV)) {  
    
    ? Note
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_m, A_alpha_m, k_alpha_m, d_alpha_m, A_beta_m, k_beta_m, d_beta_m, A_tau_m, k_tau_m, d_tau_m, A_inf_m, k_inf_m, d_inf_m
        
    TABLE minf, mtau
 DEPEND celsius
 FROM -100 TO 100 WITH 400
    
    
    UNITSOFF
    
    ? There is a Q10 factor which will alter the tau of the gates 
                 

    temp_adj_m = 3^((celsius - 24)/10)
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_m = 1
    k_alpha_m = 0.055
    d_alpha_m = -50
     
    
    alpha = A_alpha_m * exp((v - d_alpha_m) * k_alpha_m)
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_m = 1
    k_beta_m = 0.0275
    d_beta_m = -50
     
    
    beta = A_beta_m * exp((v - d_beta_m) * k_beta_m)
    
         

    ? Found a generic form of the rate equation for tau, using expression
                    tau = beta/(0.0035 *( 1 +alpha))
        
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_m = 1
    k_inf_m = -0.1
    d_inf_m = 21
     
    
    inf = A_inf_m / (exp((v - d_inf_m) * k_inf_m) + 1)
    
    minf = inf
          
       
    
    ?     *** Finished rate equations for gate
    
             

}


UNITSON