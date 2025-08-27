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
      

    SUFFIX KCa3_ChannelML_new
    USEION k READ ek WRITE ik VALENCE 1  ? reversal potential of ion is read, outgoing current is written
           
        
    USEION ca READ cai VALENCE 2 ? internal concentration of ion is read

    
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE zinf, ztau
    
}

PARAMETER { 
      

    gmax = 0.0036 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
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
    zinf
    ztau (ms)
    
}

BREAKPOINT { SOLVE states METHOD derivimplicit     

    gion = gmax*((1*m)^1)*((1*z)^1)      

    ik = gion*(v - ek)
            

}



INITIAL {
    
    ek = -80
        
    settables(v,cai)
    m = minf
           
    z = zinf
           
    
    
}
    
STATE {
    m
    z
    
}

DERIVATIVE states {
    settables(v,cai)
    m' = (minf - m)/mtau
    z' = (zinf - z)/ztau
    
}

PROCEDURE settables(v(mV), cai(mM)) {  
    
    ? Note
    LOCAL  alpha, beta, tau, inf, gamma, zeta, ca_conc, temp_adj_m, A_alpha_m, B_alpha_m, Vhalf_alpha_m, A_beta_m, B_beta_m, Vhalf_beta_m, temp_adj_z, A_alpha_z, B_alpha_z, Vhalf_alpha_z, A_beta_z, B_beta_z, Vhalf_beta_z
    
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_z = 1
    
    ? Gate depends on the concentration of ca
    ca_conc = cai ? In NEURON, the variable for the concentration  of ca is cai
    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for alpha, using expression
    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but the ChannelML file is in Physiological Units...
    ca_conc = ca_conc / 1000000
    alpha = (exp ((v-65)/27))
        
    ? Resetting concentration...
    ca_conc = ca_conc * 1000000
    
     
    ? Found a generic form of the rate equation for beta, using expression
    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but the ChannelML file is in Physiological Units...
    ca_conc = ca_conc / 1000000
    beta = 0.008
        
    ? Resetting concentration...
    ca_conc = ca_conc * 1000000
    
    mtau = 1/(temp_adj_m*(alpha + beta))
    minf = alpha/(alpha + beta)
          
       
    
    ?     *** Finished rate equations for gate
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate
         
    ? Found a generic form of the rate equation for alpha, using expression
    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but the ChannelML file is in Physiological Units...
    ca_conc = ca_conc / 1000000
    alpha = (500.0*(0.015 - (ca_conc*1e6)))/( (exp ((0.015 - (ca_conc*1e6))/0.0013)) -1)
        
    ? Resetting concentration...
    ca_conc = ca_conc * 1000000
    
     
    ? Found a generic form of the rate equation for beta, using expression
    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but the ChannelML file is in Physiological Units...
    ca_conc = ca_conc / 1000000
    beta = 0.0021
        
    ? Resetting concentration...
    ca_conc = ca_conc * 1000000
    
    ztau = 1/(temp_adj_z*(alpha + beta))
    zinf = alpha/(alpha + beta)
          
       
    
    ?     *** Finished rate equations for gate
    

         

}


UNITSON