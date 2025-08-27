INDEPENDENT {
    t FROM 0 TO 1 WITH 1 (ms)
}

NEURON {
    SUFFIX icannov
    USEION ca READ cai, ica 
    RANGE depth,taur,erev
    RANGE gbar,itrpm4, concrelease
    RANGE beta,taumin
    POINTER jip3p
    NONSPECIFIC_CURRENT itrpm4
}

UNITS {
    (mA)=(milliamp)
    (mV)=(millivolt)
    (molar)=(1/liter)
    (mM)=(millimolar)
    (um)=(micron)
    (msM)=(ms mM)
    FARADAY=(faraday) (coulomb)
	PI      = (pi)       (1)
}

PARAMETER {
    v (mV)
    depth=0.0125 (um)         
    taur=80 (ms)           
    erev=0 (mV)             
    cai (mM)                
	ica       (mA/cm2)
    gbar=0.0001 (mho/cm2)
    
    taumin=0.1 (ms)         
    concrelease=500
    Kd = 87e-3 (mM)		
}

STATE {
    can (mM) 
    Po
}

ASSIGNED {
	jip3p (mM/ms)
    ican (mA/cm2)
    drive_channel (mM/ms)
    itrpm4 (mA/cm2)
    Po_inf
    Tau (ms)
	diam      (um)
    
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    itrpm4=gbar*Po*(v-erev)
}

DERIVATIVE states {
    evaluate_fct(v,can)
    Po'=(Po_inf-Po)/(Tau)
	drive_channel =  - (concrelease) * ica / (2 * FARADAY * depth)
    if (drive_channel<=0.0) {drive_channel=0.0}             
    can'=drive_channel+(cai-can)/(taur)
}

FUNCTION MyExp(x) {
    if (x<-50) {MyExp=0}
    else if (x>50) {MyExp=exp(50)}
    else {MyExp=exp(x)}
}

UNITSOFF

INITIAL {
    
    
    can=cai
    evaluate_fct(v,can)
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {
    LOCAL alpha, alpha2, beta
    
    
    
    alpha=0.0057*MyExp(0.0060*-60)
    beta=0.033*MyExp(-0.019*-60)

    alpha2=alpha/(1+(Kd/cai))
    Po_inf=alpha2/(alpha2+beta)
    Tau=1/(alpha2+beta)
    if (Tau<taumin) {Tau=taumin}                        
}

UNITSON