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
      

    SUFFIX BK_IAMC_ChannelML
    USEION k READ ek WRITE ik VALENCE 1  ? reversal potential of ion is read, outgoing current is written
           
        
    USEION ca READ cai VALENCE 2 ? internal concentration of ion is read

    
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE ninf, ntau
    
}

PARAMETER { 
      

    gmax = 0.000014999999999999999 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
          

    ? Reversal potential of k
    ek (mV)
    ? The outward flow of ion
    ik (mA/cm2)
          

    ? The internal concentration of ion
    cai (mM)   
    
    
    gion (S/cm2)
    minf
    mtau (ms)
    ninf
    ntau (ms)
    
}

BREAKPOINT { 
    SOLVE states METHOD derivimplicit
    gion = gmax * (m
^1) * (n
^1)      

    ik = gion*(v - ek)
            

}



INITIAL {
    
    ek = -86.5
        
    settables(v,cai)
    m = minf
        n = ninf
        
    
}
    
STATE {
    m
    n
    
}



DERIVATIVE states {
    settables(v,cai)
    m' = (minf - m)/mtau
            n' = (ninf - n)/ntau
            

}

PROCEDURE settables(v(mV), cai(mM)) {  
    
    ? Note
    LOCAL  alpha, beta, tau, inf, gamma, zeta, ca_conc
, temp_adj_m
, temp_adj_n
    
    UNITSOFF
    
    ? There is a Q10 factor which will alter the tau of the gates 
                 

    temp_adj_m = 3^((celsius - 17.350264793)/10)     

    temp_adj_n = 3^((celsius - 17.350264793)/10)
    
    ? There is a voltage offset of 0.010. This will shift the dependency of the rate equations 
    v = v - (10)
    
    ? Gate depends on the concentration of ca
    ca_conc = cai ? In NEURON, the variable for the concentration  of ca is cai
    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for alpha, using expression
    
    ? Note
    
    v = v * 0.001   ? temporarily set v to units of equation...
            
    alpha = 2500/(1 + ( (1.5e-3 *(exp (-85*v))) / ca_conc))
        
    ? Set correct units of alpha for NEURON
    alpha = alpha * 0.001 
    
    v = v * 1000   ? reset v
        
     
    ? Found a generic form of the rate equation for beta, using expression
    
    ? Note
    
    v = v * 0.001   ? temporarily set v to units of equation...
            
    beta = 1500/(1 + (ca_conc / (1.5e-4 * (exp (-77*v)))))
        
    ? Set correct units of beta for NEURON
    beta = beta * 0.001 
    
    v = v * 1000   ? reset v
        
    mtau = 1/(temp_adj_m*(alpha + beta))
    minf = alpha/(alpha + beta)
    


    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for tau, using expression
    
    ? Note
    
    v = v * 0.001   ? temporarily set v to units of equation...
            
    tau = 0.005
        
    ? Set correct units of tau for NEURON
    tau = tau * 1000 
    
    v = v * 1000   ? reset v
        
    ntau = tau/temp_adj_n
     
    ? Found a generic form of the rate equation for inf, using expression
    
    ? Note
    
    v = v * 0.001   ? temporarily set v to units of equation...
            
    inf = 1
         
    
    v = v * 1000   ? reset v
        
    ninf = inf
    


    ?     *** Finished rate equations for gate
    

         

}


UNITSON