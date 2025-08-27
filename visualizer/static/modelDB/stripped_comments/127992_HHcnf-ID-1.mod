UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
} 


NEURON {
    SUFFIX HHcnf
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE el, gl
    RANGE gnabar, gkbar
    RANGE gna, gk
    RANGE gamma_na, gamma_k
    RANGE m_inf, h_inf, n_inf
    RANGE tau_m, tau_h, tau_n, tau_y, tau_z
    RANGE var_y, var_z
    RANGE noise_y, noise_z
    RANGE Nna, Nk
    RANGE seed    
    
    RANGE mu_y, mu_z
    THREADSAFE
} 


PARAMETER {
    gnabar  = 0.12   (S/cm2)
    gkbar   = 0.036  (S/cm2)
    gl      = 0.0003 (S/cm2)  
    el       = -54.3 (mV)       
    
    gamma_na = 10  (pS)		
    gamma_k  = 10  (pS)		
    seed = 5061983              
} 


STATE {
    m h n yy zz
} 


ASSIGNED {
    ina        (mA/cm2)
    ik         (mA/cm2)
    il         (mA/cm2)
    gna        (S/cm2)
    gk         (S/cm2)
    ena        (mV)
    ek         (mV)
    
    dt         (ms)
    area       (um2)
    celsius    (degC)
    v          (mV)
        
    Nna  (1) 
    Nk   (1) 
    
    m_inf h_inf n_inf
    noise_y noise_z
    var_y (ms2) var_z (ms2)
    tau_m (ms) tau_h (ms) tau_n (ms) tau_y (ms) tau_z (ms)
    
    mu_y mu_z
} 

INITIAL {
    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   
    
    rates(v)
    m = m_inf
    h = h_inf
    n = n_inf
    yy = 0.
    zz = 0.
    set_seed(seed)
} 


BREAKPOINT {
    SOLVE states
    gna = gnabar * (m*m*m*h + zz)
    if (gna < 0) {
	gna = 0
    }
    gk = gkbar * (n*n*n*n + yy)
    if (gk < 0) {
	gk = 0
    }
    ina = gna * (v - ena)
    ik  = gk * (v - ek)
    il  = gl * (v - el)
} 


PROCEDURE states() {
    rates(v)
    m = m + dt * (m_inf-m)/tau_m
    h = h + dt * (h_inf-h)/tau_h
    n = n + dt * (n_inf-n)/tau_n
    
    yy = yy*mu_y + noise_y
    zz = zz*mu_z + noise_z
    
    
    
    
    
    VERBATIM
    return 0;
    ENDVERBATIM
} 


PROCEDURE rates(Vm (mV)) { 
    LOCAL a,b,m3_inf,n4_inf,q10,sum,var_z1,var_z2,var_z3,var_z4,var_z5,var_z6,var_z7,tau_z1,tau_z2,tau_z3,tau_z4,tau_z5,tau_z6,tau_z7,var_y1,var_y2,var_y3,var_y4,one_minus_m,one_minus_h,one_minus_n
    
    UNITSOFF
    
    q10 = 3.^((celsius-6.3)/10.)
    
    
    a = alpham(Vm)
    b = betam(Vm)
    sum = a+b
    tau_m = 1. / (q10*sum)
    m_inf = a / sum
    one_minus_m = 1. - m_inf
    m3_inf = m_inf*m_inf*m_inf

    
    a = alphah(Vm)
    b = betah(Vm)
    sum = a+b
    tau_h = 1. / (q10*sum)
    h_inf = a / sum
    one_minus_h = 1. - h_inf
    
    tau_z1 = tau_h
    tau_z2 = tau_m
    tau_z3 = tau_m/2
    tau_z4 = tau_m/3
    tau_z5 = tau_m*tau_h/(tau_m+tau_h)
    tau_z6 = tau_m*tau_h/(tau_m+2*tau_h)
    tau_z7 = tau_m*tau_h/(tau_m+3*tau_h)
    var_z1 = 1.0 / Nna * m3_inf*m3_inf*h_inf * one_minus_h
    var_z2 = 3.0 / Nna * m3_inf*m_inf*m_inf*h_inf*h_inf * one_minus_m
    var_z3 = 3.0 / Nna * m3_inf*m_inf*h_inf*h_inf * one_minus_m*one_minus_m
    var_z4 = 1.0 / Nna * m3_inf*h_inf*h_inf * one_minus_m*one_minus_m*one_minus_m
    var_z5 = 3.0 / Nna * m3_inf*m_inf*m_inf*h_inf * one_minus_m*one_minus_h
    var_z6 = 3.0 / Nna * m3_inf*m_inf*h_inf * one_minus_m*one_minus_m*one_minus_h
    var_z7 = 1.0 / Nna * m3_inf*h_inf * one_minus_m*one_minus_m*one_minus_m*one_minus_h
    
    tau_z = (var_z1+var_z2+var_z3+var_z4+var_z5+var_z6+var_z7) / (var_z1/tau_z1+var_z2/tau_z2+var_z3/tau_z3+var_z4/tau_z4+var_z5/tau_z5+var_z6/tau_z6+var_z7/tau_z7)
    var_z = 1.0 / Nna * m3_inf * h_inf * (1.0 - m3_inf*h_inf)
     
    
    mu_z = exp(-dt/tau_z)
    noise_z = sqrt(var_z * (1-mu_z*mu_z)) * normrand(0,1)
    
    
    
    
    a = alphan(Vm)
    b = betan(Vm)
    sum = a+b
    tau_n = 1. / (q10*sum)
    n_inf = a / sum
    one_minus_n = 1. - n_inf
    n4_inf = n_inf * n_inf * n_inf * n_inf
    
    var_y1 = 4.0/Nk * n4_inf*n_inf*n_inf*n_inf * one_minus_n
    var_y2 = 6.0/Nk * n4_inf*n_inf*n_inf * one_minus_n*one_minus_n
    var_y3 = 4.0/Nk * n4_inf*n_inf * one_minus_n*one_minus_n*one_minus_n
    var_y4 = 1.0/Nk * n4_inf * one_minus_n*one_minus_n*one_minus_n*one_minus_n

    tau_y = (var_y1+var_y2+var_y3+var_y4) / (var_y1/tau_n + var_y2/tau_n/2 + var_y3/tau_n/3 + var_y4/tau_n/4)
    var_y = 1.0/Nk * n4_inf * (1.0-n4_inf)
    
    
    mu_y = exp(-dt/tau_y)
    noise_y = sqrt(var_y * (1-mu_y*mu_y)) * normrand(0,1)
    
    
    
    UNITSON
}

PROCEDURE consistency(){
    if (m > 1) {
	m = 1
    }
    if (m < 0) {
	m = 0
    }
    if (h > 1) {
        h = 1
    }
    if (h < 0) {
	h = 0
    }
    if (n > 1) {
	n = 1
    }
    if (n < 0) {
	n = 0
    } 
}


FUNCTION vtrap(x,y) {  
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}


FUNCTION alpham(Vm (mV)) (/ms) {
    UNITSOFF
    alpham = .1 * vtrap(-(Vm+40),10)
    UNITSON
}


FUNCTION betam(Vm (mV)) (/ms) {
    UNITSOFF
    betam =  4 * exp(-(Vm+65)/18)
    UNITSON
}


FUNCTION alphah(Vm (mV)) (/ms) {
    UNITSOFF
    alphah = .07 * exp(-(Vm+65)/20)
    UNITSON
}


FUNCTION betah(Vm (mV)) (/ms) {
    UNITSOFF
    betah = 1 / (exp(-(Vm+35)/10) + 1)
    UNITSON
}


FUNCTION alphan(Vm (mV)) (/ms) {
    UNITSOFF
    alphan = .01*vtrap(-(Vm+55),10) 
    UNITSON
}


FUNCTION betan(Vm (mV)) (/ms) {
    UNITSOFF
    betan = .125*exp(-(Vm+65)/80)
    UNITSON
}