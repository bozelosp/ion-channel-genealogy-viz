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
      

    SUFFIX LCa3_mit_usb_ChannelML
    USEION ca READ eca WRITE ica VALENCE 2 ? reversal potential of ion is read, outgoing current is written
            
    RANGE gmax, gion
    
    RANGE minf, mtau
    RANGE hinf, htau
}

PARAMETER { 
      

    gmax = 0.012 (S/cm2) ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
    
    ? Reversal potential of ca
    eca (mV)
    ? The outward flow of ion
    ica (mA/cm2)
            
    
    gion (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
         

    gion = gmax*((1*m)^1)*((1*h)^1)
    ica = gion*(v - eca)
                

}



INITIAL {
    eca = 70
        
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
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_m, A_alpha_m, k_alpha_m, d_alpha_m, A_beta_m, k_beta_m, d_beta_m, temp_adj_h, A_alpha_h, k_alpha_h, d_alpha_h, A_beta_h, k_beta_h, d_beta_h
        
    TABLE minf, mtau,hinf, htau
 DEPEND celsius
 FROM -100 TO 50 WITH 3000
    
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_m = 7500
    k_alpha_m = -142.85714285714286
    d_alpha_m = 0.013
    
    ? Unit system in ChannelML file is SI units, therefore need to 
    ? convert these to NEURON quanities...
                        A_alpha_m = A_alpha_m * 0.0010   ? 1/ms
    k_alpha_m = k_alpha_m * 0.0010   ? mV
    d_alpha_m = d_alpha_m * 1000   ? mV
          
                     
    
    alpha = A_alpha_m / (exp((v - d_alpha_m) * k_alpha_m) + 1)
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_m = 1650
    k_beta_m = 250
    d_beta_m = 0.014
    
    ? Unit system in ChannelML file is SI units, therefore need to 
    ? convert these to NEURON quanities...
                        A_beta_m = A_beta_m * 0.0010   ? 1/ms
    k_beta_m = k_beta_m * 0.0010   ? mV
    d_beta_m = d_beta_m * 1000   ? mV
          
                     
    
    beta = A_beta_m / (exp((v - d_beta_m) * k_beta_m) + 1)
    
    mtau = 1/(temp_adj_m*(alpha + beta))
    minf = alpha/(alpha + beta)
          
       
    
    ?     *** Finished rate equations for gate
    
        
        
    ?      ***  Adding rate equations for gate
        
    ? Found a parameterised form of rate equation for alpha, using expression
    A_alpha_h = 6.800
    k_alpha_h = 83.333333333333
    d_alpha_h = -0.030
    
    ? Unit system in ChannelML file is SI units, therefore need to 
    ? convert these to NEURON quanities...
                        A_alpha_h = A_alpha_h * 0.0010   ? 1/ms
    k_alpha_h = k_alpha_h * 0.0010   ? mV
    d_alpha_h = d_alpha_h * 1000   ? mV
          
                     
    
    alpha = A_alpha_h / (exp((v - d_alpha_h) * k_alpha_h) + 1)
    
    
    ? Found a parameterised form of rate equation for beta, using expression
    A_beta_h = 60
    k_beta_h = -90.90909090909
    d_beta_h = 0.0
    
    ? Unit system in ChannelML file is SI units, therefore need to 
    ? convert these to NEURON quanities...
                        A_beta_h = A_beta_h * 0.0010   ? 1/ms
    k_beta_h = k_beta_h * 0.0010   ? mV
    d_beta_h = d_beta_h * 1000   ? mV
          
                     
    
    beta = A_beta_h / (exp((v - d_beta_h) * k_beta_h) + 1)
    
    htau = 1/(temp_adj_h*(alpha + beta))
    hinf = alpha/(alpha + beta)
          
       
    
    ?     *** Finished rate equations for gate
    
             

}


UNITSON