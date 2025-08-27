NEURON {
    SUFFIX iM
    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!
    
    RANGE gion                           
    RANGE gmax                              
    RANGE conductance                       
    
    RANGE g                                 
    
    RANGE fopen                             
    RANGE m_instances                       
    
    RANGE m_tau                             
    GLOBAL m_min    
    RANGE m_inf                             
    
    RANGE m_rateScale                       
    
    RANGE m_fcond                           
    RANGE m_timeCourse_TIME_SCALE           
    RANGE m_timeCourse_VOLT_SCALE           
    
    RANGE m_timeCourse_t                    
    RANGE m_steadyState_rate                
    RANGE m_steadyState_midpoint            
    RANGE m_steadyState_scale               
    
    RANGE m_steadyState_x                   
    RANGE m_q10Settings_q10Factor           
    RANGE m_q10Settings_experimentalTemp    
    RANGE m_q10Settings_TENDEGREES          
    
    RANGE m_q10Settings_q10                 
    RANGE m_timeCourse_V                    
    RANGE m_tauUnscaled                     
    RANGE conductanceScale                  
    RANGE fopen0                            
    RANGE gmin
    RANGE shift
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
    gmin = 0  (S/cm2)
    conductance = 1.0E-5 (uS)
    m_instances = 1 
    m_timeCourse_TIME_SCALE = 1 (ms)
    m_timeCourse_VOLT_SCALE = 1 (mV)
    m_steadyState_rate = 1 
    m_steadyState_midpoint = -36.7 (mV)
    m_steadyState_scale = 9.48 (mV)
    m_q10Settings_q10Factor = 2.6
    m_q10Settings_experimentalTemp = 308.15 (K)
    m_q10Settings_TENDEGREES = 10 (K)
    shift=3 (mV)
    m_min=0
}

ASSIGNED {
    
    gion   (S/cm2)                          
    v (mV)
    celsius (degC)
    temperature (K)
    ek (mV)
    ik (mA/cm2)
    
    
    m_timeCourse_V                         
    
    m_timeCourse_t (ms)                    
    
    m_steadyState_x                        
    
    m_q10Settings_q10                      
    
    m_rateScale                            
    
    m_fcond                                
    
    m_inf                                  
    
    m_tauUnscaled (ms)                     
    
    m_tau (ms)                             
    
    conductanceScale                       
    
    fopen0                                 
    
    fopen                                  
    
    g (uS)                                 
    rate_m_q (/ms)
    
}

STATE {
    m_q  
    
}

INITIAL {
    ek = -90.0
    
    temperature = celsius + 273.15
    
    rates()
    rates() ? To ensure correct initialisation.
    
    m_q = m_inf
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    
    ? DerivedVariable is based on path
    ? Path not present in component, using factor
    
    conductanceScale = 1 
    
    ? DerivedVariable is based on path
    ? multiply applied to all instances of fcond in
    fopen0 = m_fcond ? path based, prefix = 
    
    fopen = conductanceScale  *  fopen0 ? evaluable
    g = conductance  *  fopen   ? evaluable
    gion = gmax * (m_min+fopen) + gmin 
    
    ik = gion * (v - ek)
    
}

DERIVATIVE states {
    rates()
    m_q' = rate_m_q 
    
}

PROCEDURE rates() {
    
    m_timeCourse_V = v /  m_timeCourse_VOLT_SCALE ? evaluable
    m_timeCourse_t = (13.4+26.3*exp(    -( ( m_timeCourse_V -29.7+shift)/30.3 ) * ( ( m_timeCourse_V -29.7+shift)/30.3 )   )) *  m_timeCourse_TIME_SCALE ? evaluable
    m_steadyState_x = m_steadyState_rate  / (1 + exp(0 - (v -  m_steadyState_midpoint +shift)/ m_steadyState_scale )) ? evaluable
    m_q10Settings_q10 = m_q10Settings_q10Factor ^((temperature -  m_q10Settings_experimentalTemp )/ m_q10Settings_TENDEGREES ) ? evaluable
    ? DerivedVariable is based on path
    ? multiply applied to all instances of q10 in
    m_rateScale = m_q10Settings_q10 ? path based, prefix = m_
    
    m_fcond = m_q ^ m_instances ? evaluable
    ? DerivedVariable is based on path
    m_inf = m_steadyState_x ? path based, prefix = m_
    
    ? DerivedVariable is based on path
    m_tauUnscaled = m_timeCourse_t ? path based, prefix = m_
    
    m_tau = m_tauUnscaled  /  m_rateScale ? evaluable
    
     
    rate_m_q = ( m_inf  -  m_q ) /  m_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    
     
    
}