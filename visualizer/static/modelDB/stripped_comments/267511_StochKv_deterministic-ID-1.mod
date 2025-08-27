NEURON {
    SUFFIX StochKv_deterministic
    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!
    
    RANGE gion                           
    RANGE gmax                              
    RANGE conductance                       
    
    RANGE g                                 
    
    RANGE fopen                             
    RANGE q10ConductanceScaling_q10Factor   
    RANGE q10ConductanceScaling_experimentalTemp
    
    RANGE q10ConductanceScaling_factor      
    RANGE n_instances                       
    
    RANGE n_alpha                           
    
    RANGE n_beta                            
    
    RANGE n_tau                             
    
    RANGE n_inf                             
    
    RANGE n_rateScale                       
    
    RANGE n_fcond                           
    RANGE n_reverseRate_rate                
    RANGE n_reverseRate_midpoint            
    RANGE n_reverseRate_scale               
    
    RANGE n_reverseRate_r                   
    RANGE n_forwardRate_rate                
    RANGE n_forwardRate_midpoint            
    RANGE n_forwardRate_scale               
    
    RANGE n_forwardRate_r                   
    RANGE n_q10Settings_q10Factor           
    RANGE n_q10Settings_experimentalTemp    
    RANGE n_q10Settings_TENDEGREES          
    
    RANGE n_q10Settings_q10                 
    RANGE n_reverseRate_x                   
    RANGE n_forwardRate_x                   
    RANGE conductanceScale                  
    RANGE fopenHHrates                      
    RANGE fopenHHtauInf                     
    RANGE fopenHHratesTau                   
    RANGE fopenHHratesInf                   
    RANGE fopenHHratesTauInf                
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    
}

PARAMETER {
    
    gmax = 0  (S/cm2)                       
    
    conductance = 1.0E-5 (uS)
    q10ConductanceScaling_q10Factor = 2.3 
    q10ConductanceScaling_experimentalTemp = 296.15 (K)
    n_instances = 1 
    n_reverseRate_rate = 0.018000001 (kHz)
    n_reverseRate_midpoint = -40 (mV)
    n_reverseRate_scale = -9 (mV)
    n_forwardRate_rate = 0.18 (kHz)
    n_forwardRate_midpoint = -40 (mV)
    n_forwardRate_scale = 9 (mV)
    n_q10Settings_q10Factor = 2.3 
    n_q10Settings_experimentalTemp = 296.15 (K)
    n_q10Settings_TENDEGREES = 10 (K)
}

ASSIGNED {
    
    gion   (S/cm2)                          
    v (mV)
    celsius (degC)
    temperature (K)
    ek (mV)
    ik (mA/cm2)
    
    
    q10ConductanceScaling_factor           
    
    n_reverseRate_x                        
    
    n_reverseRate_r (kHz)                  
    
    n_forwardRate_x                        
    
    n_forwardRate_r (kHz)                  
    
    n_q10Settings_q10                      
    
    n_rateScale                            
    
    n_alpha (kHz)                          
    
    n_beta (kHz)                           
    
    n_fcond                                
    
    n_inf                                  
    
    n_tau (ms)                             
    
    conductanceScale                       
    
    fopenHHrates                           
    
    fopenHHtauInf                          
    
    fopenHHratesTau                        
    
    fopenHHratesInf                        
    
    fopenHHratesTauInf                     
    
    fopen                                  
    
    g (uS)                                 
    rate_n_q (/ms)
    
}

STATE {
    n_q 
    
}

INITIAL {
    ek = -85.0
    
    temperature = celsius + 273.15
    
    rates()
    rates() ? To ensure correct initialisation.
    
    n_q = n_inf
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    
    ? DerivedVariable is based on path
    conductanceScale = q10ConductanceScaling_factor ? multiply applied to all instances of factor in
    
    ? DerivedVariable is based on path
    fopenHHrates = n_fcond ? multiply applied to all instances of fcond in
    
    ? DerivedVariable is based on path
    ? Path not present in component, using factor
    
    fopenHHtauInf = 1 
    
    ? DerivedVariable is based on path
    ? Path not present in component, using factor
    
    fopenHHratesTau = 1 
    
    ? DerivedVariable is based on path
    ? Path not present in component, using factor
    
    fopenHHratesInf = 1 
    
    ? DerivedVariable is based on path
    ? Path not present in component, using factor
    
    fopenHHratesTauInf = 1 
    
    fopen = conductanceScale  *  fopenHHrates  *  fopenHHtauInf  *  fopenHHratesTau  *  fopenHHratesInf  *  fopenHHratesTauInf ? evaluable
    g = conductance  *  fopen ? evaluable
    gion = gmax * fopen 
    
    ik = gion * (v - ek)
    
}

DERIVATIVE states {
    rates()
    n_q' = rate_n_q 
    
}

PROCEDURE rates() {
    
    q10ConductanceScaling_factor = q10ConductanceScaling_q10Factor ^((temperature -  q10ConductanceScaling_experimentalTemp )/10) ? evaluable
    n_reverseRate_x = (v -  n_reverseRate_midpoint ) /  n_reverseRate_scale ? evaluable
    if (n_reverseRate_x  != 0)  { 
        n_reverseRate_r = n_reverseRate_rate  *  n_reverseRate_x  / (1 - exp(0 -  n_reverseRate_x )) ? evaluable cdv
    } else if (n_reverseRate_x  == 0)  { 
        n_reverseRate_r = n_reverseRate_rate ? evaluable cdv
    }
    
    n_forwardRate_x = (v -  n_forwardRate_midpoint ) /  n_forwardRate_scale ? evaluable
    if (n_forwardRate_x  != 0)  { 
        n_forwardRate_r = n_forwardRate_rate  *  n_forwardRate_x  / (1 - exp(0 -  n_forwardRate_x )) ? evaluable cdv
    } else if (n_forwardRate_x  == 0)  { 
        n_forwardRate_r = n_forwardRate_rate ? evaluable cdv
    }
    
    n_q10Settings_q10 = n_q10Settings_q10Factor ^((temperature -  n_q10Settings_experimentalTemp )/ n_q10Settings_TENDEGREES ) ? evaluable
    ? DerivedVariable is based on path
    n_rateScale = n_q10Settings_q10 ? multiply applied to all instances of q10 in
    
    ? DerivedVariable is based on path
    n_alpha = n_forwardRate_r ? path based
    
    ? DerivedVariable is based on path
    n_beta = n_reverseRate_r ? path based
    
    n_fcond = n_q ^ n_instances ? evaluable
    n_inf = n_alpha /( n_alpha + n_beta ) ? evaluable
    n_tau = 1/(( n_alpha + n_beta ) *  n_rateScale ) ? evaluable
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    rate_n_q = ( n_inf  -  n_q ) /  n_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    
     
    
}