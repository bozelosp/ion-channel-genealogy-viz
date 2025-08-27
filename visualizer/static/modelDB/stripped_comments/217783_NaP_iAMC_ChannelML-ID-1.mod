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
      

    SUFFIX NaP_iAMC_ChannelML
    USEION na READ ena WRITE ina VALENCE 1 ? reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE hinf, htau
    
    RANGE ninf, ntau
    
}

PARAMETER { 
      

    gmax = 0.000059999999999999995 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
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
    ninf
    ntau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * ((1*m)
^3) * ((1*h)
^1) * ((1*n)
^1)      

    ina = gion*(v - ena)
            

}



INITIAL {
    
    ena = 67
        
    rates(v)
    m = minf
        h = hinf
        n = ninf
        
    
}
    
STATE {
    m
    h
    n
    
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
            h' = (hinf - h)/htau
            n' = (ninf - n)/ntau
            

}

PROCEDURE rates(v(mV)) {  
    
    ? Note
    LOCAL  alpha, beta, tau, inf, gamma, zeta
, temp_adj_m,
         A_inf_m, B_inf_m, Vhalf_inf_m
, temp_adj_h,
         A_inf_h, B_inf_h, Vhalf_inf_h
, temp_adj_n,
         A_inf_n, B_inf_n, Vhalf_inf_n
    
    TABLE minf, mtau,hinf, htau,ninf, ntau
 DEPEND celsius FROM -100 TO 100 WITH 2000
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    temp_adj_n = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    tau = (1+(4 * (exp(0 - ((v + 50)/20)^2))))
        
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_m = 0.499622025796
    B_inf_m = -4.9
    Vhalf_inf_m = -59.0 
    inf = A_inf_m / (exp((v - Vhalf_inf_m) / B_inf_m) + 1)
    
    minf = inf
    


    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    tau = (5000+(16000 * (exp(0 - ((v + 50)/20)^2))))
        
    htau = tau/temp_adj_h
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_h = 0.499622025796
    B_inf_h = 4.9
    Vhalf_inf_h = -59.0 
    inf = A_inf_h / (exp((v - Vhalf_inf_h) / B_inf_h) + 1)
    
    hinf = inf
    


    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    tau = 1
        
    ntau = tau/temp_adj_n
    
    ? Found a parameterised form of rate equation for inf, using expression
    A_inf_n = 0.499622025796
    B_inf_n = 4.9
    Vhalf_inf_n = -59.0 
    inf = A_inf_n / (exp((v - Vhalf_inf_n) / B_inf_n) + 1)
    
    ninf = inf
    


    ?     *** Finished rate equations for gate
    

         

}


UNITSON