INDEPENDENT {
    t FROM 0 TO 1 WITH 1 (ms)
}

NEURON {
    SUFFIX ican
    USEION ca READ cai 
    USEION na WRITE ina
    RANGE depth,taur,frac,ica,ican
    RANGE gbar,m_inf,tau_m,in,mystart
    RANGE beta,cac,taumin
    POINTER cip3p
}

UNITS {
    (mA)=(milliamp)
    (mV)=(millivolt)
    (molar)=(1/liter)
    (mM)=(millimolar)
    (um)=(micron)
    (msM)=(ms mM)
    FARADAY=(faraday) (coulomb)
}

PARAMETER {
    v (mV)
    depth=0.08 (um)         
    frac=0
    taur=100 (ms)           
    
    en=-20 (mV)             
    eca=140 (mV)            
    cai (mM)                
    gbar=0.0001 (mho/cm2)
    beta=0.0001 (1/ms)      
    cac=0.15 (mM)	
    
    taumin=0.1 (ms)         
}

STATE {
    can (mM) 
    m
}

ASSIGNED {
	cip3p (mA/cm2)
    ica (mA/cm2)
    ican (mA/cm2)
    drive_channel (mM/ms)
    in (mA/cm2)
    ina (mA/cm2)
    m_inf
    tau_m (ms)
    tadj
    
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    in=gbar*m*m*(v-en)
    ican=gbar*frac*m*m*(v-eca)
    ina=0.7*in
}

DERIVATIVE states {
    evaluate_fct(v,can)
    m'=(m_inf-m)/tau_m
    drive_channel=-(10000)*(ican-100*cip3p)/(2*FARADAY*depth) 
    if (drive_channel<=0.0) {drive_channel=0.0}             
    can'=drive_channel/18+(cai-can)/(taur)
}

UNITSOFF

INITIAL {
    
    
    can=cai
    tadj=3.0^((celsius-22.0)/10)
    evaluate_fct(v,can)
    m=m_inf
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {
    LOCAL alpha2
    alpha2=beta*(cai/cac)^2
    if (cai<1.09e-4) {alpha2=0}
    tau_m=1/(alpha2+beta)/tadj
    m_inf=alpha2/(alpha2+beta)
    if (tau_m<taumin) {tau_m=taumin}                        
}

UNITSON