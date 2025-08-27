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
      

    SUFFIX KA_iAMC_ChannelML
    USEION k READ ek WRITE ik VALENCE 1  ? reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE hinf, htau
    
}

PARAMETER { 
      

    gmax = 0.004 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
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
        
    gion = gmax * (m
^1) * (h
^1)      

    ik = gion*(v - ek)
            

}



INITIAL {
    
    ek = -86.5
        
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
, temp_adj_m
, temp_adj_h,
         A_inf_h, B_inf_h, Vhalf_inf_h
    
    TABLE minf, mtau,hinf, htau
 DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    tau = (1+4*exp(-((v-32)/50)^2))
        
    mtau = tau/temp_adj_m
     
    ? Found a generic form of the rate equation for inf, using expression
    
    
    if (v < -50 ) {
        inf =  0 
    } else {
        inf =  1 / (1 + exp(0 - (v + 25.7)/4.4))
    }
    minf = inf
    


    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    tau = (1+100*exp(-((v-10)/40)^2))
        
    htau = tau/temp_adj_h
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_h = 1
    B_inf_h = 4.4
    Vhalf_inf_h = -25 
    inf = A_inf_h / (exp((v - Vhalf_inf_h) / B_inf_h) + 1)
    
    hinf = inf
    


    ?     *** Finished rate equations for gate
    

         

}


UNITSON