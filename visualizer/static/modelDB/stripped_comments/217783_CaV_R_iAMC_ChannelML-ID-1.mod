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
      

    SUFFIX CaV_R_iAMC_ChannelML
    USEION ca WRITE ica VALENCE 2 ?  outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE hinf, htau
    
}

PARAMETER { 
      

    gmax = 0.00015 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
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
        
    gion = gmax * (m
^2) * (h
^1)      

    ica = gion*(v - eca)
            

}



INITIAL {
    
    eca = 80
        
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
    LOCAL  alpha, beta, tau, inf, gamma, zeta
, temp_adj_m,
         A_inf_m, B_inf_m, Vhalf_inf_m
, temp_adj_h
    
    TABLE minf, mtau,hinf, htau
 DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    
    
    if (v < -30 ) {
        tau =  28.4118 
    } else {
        tau =  3.1738 + (25.238 * (exp(-1 * ((v + 30)/17.498))))
    }
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_m = 1
    B_inf_m = -2.0914
    Vhalf_inf_m = -38.037 
    inf = A_inf_m / (exp((v - Vhalf_inf_m) / B_inf_m) + 1)
    
    minf = inf
    


    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    
    
    if (v < -30 ) {
        tau =  21.0638148543 
    } else {
        tau =  10.8 + (3.0 * (exp(-1 * ((v+20)/8.13))))
    }
    htau = tau/temp_adj_h
     
    ? Found a generic form of the rate equation for inf, using expression
    inf = ((1/(1+(exp(-1 * (v-(-38.037))/-2.0914)))) + (0.6928/(1+(exp(-1 * (v-(-38.037))/2.0914)))))
        
    hinf = inf
    


    ?     *** Finished rate equations for gate
    

         

}


UNITSON