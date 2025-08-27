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
      

    SUFFIX Kdr_iAMC_Fig10Hii_ChannelML
    USEION k READ ek WRITE ik VALENCE 1 ? reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE ninf, ntau
    
}

PARAMETER { 
      

    gmax = 0.001 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
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
    ninf
    ntau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * ((1*m)
^1) * (n
^1)      

    ik = gion*(v - ek)
            

}



INITIAL {
    
    ek = -86.5
        
    rates(v)
    m = minf
        n = ninf
        
    
}
    
STATE {
    m
    n
    
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
            n' = (ninf - n)/ntau
            

}

PROCEDURE rates(v(mV)) {  
    
    ? Note
    LOCAL  alpha, beta, tau, inf, gamma, zeta
, temp_adj_m
, temp_adj_n
    
    TABLE minf, mtau,ninf, ntau
 DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_n = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    
    
    if (v < -66.67 ) {
        tau =  230.256604897 
    } else {
        tau =  29.156*exp(-(v+55)/5.8842) + 18.394
    }
    mtau = tau/temp_adj_m
     
    ? Found a generic form of the rate equation for inf, using expression
    
    
    if (v < -66.67 ) {
        inf =  0 
    } else {
        inf =  0.005*v + 0.35
    }
    minf = inf
    


    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    tau = 4
        
    ntau = tau/temp_adj_n
     
    ? Found a generic form of the rate equation for inf, using expression
    
    
    if (v < 0 ) {
        inf =  0 
    } else {
        inf =  ((-0.006*v) + 0.40)
    }
    ninf = inf
    


    ?     *** Finished rate equations for gate
    

         

}


UNITSON