NEURON {
    POINT_PROCESS Wghkampa_preML
    USEION na WRITE ina
    USEION k WRITE ik
    USEION ca READ cai	
    USEION glut READ gluti WRITE iglut,gluti VALENCE 0
    
    RANGE taur, taud
    RANGE iampa,winit, iglut
    RANGE P, Pmax, lr
    
    
    RANGE g_factor, I, U_SE_factor, glut_factor
    RANGE U_SE, U_SE_init, tau_takeover, tau_in
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {
    taur=2 		(ms) <1e-9,1e9>
    taud = 10 	(ms) <1e-9,1e9>
    nai = 18	(mM)	
    nao = 140	(mM)
    ki = 140	(mM)	
    ko = 5		(mM)
    cai			(mM)
    celsius		(degC)
    Pmax=1e-6   (cm/s)	
    alpha1=0.35	
    beta1=80
    alpha2=0.55
    beta2=80
    
    
    
    tau_rec = 0.8 (ms)
    tau_in = 3 (ms)
    U_SE_init = 0.1 (1)
    U_SE_factor = 0 (1)
    tau_takeover = 0.1 (ms)
    
    
    glut_factor = 40 (1) 
    g_factor = 1 (1)
}


ASSIGNED {
    ina     (nA)
    ik      (nA)
    v (mV)
    P (cm/s)
    factor
    iampa	(nA)
    lr
    Area (cm2)
    U_SE (1)
    
    I (1)
    
    
    iglut (nA)
}


STATE {
    A (cm/s)
    B (cm/s)
    
    x  (1)
    y_rel  (1)
    z  (1)    
    gluti (mM)
}

INITIAL {
    LOCAL tp
    if (taur/taud > .9999) {
	taur = .9999*taud
    }
    A = 0
    B = 0
    tp = (taur*taud)/(taud - taur) * log(taud/taur)
    factor = -exp(-tp/taur) + exp(-tp/taud)
    factor = 1/factor
    Area=1
    U_SE = U_SE_init
    
    
    
    x = 1 (1)
    y_rel = 0 (1)
    z = 0 (1)
    I = 0 (1)

    gluti = 0 (mM)
}


BREAKPOINT {
    SOLVE state METHOD cnexp
    P=B-A
    
    
    
    
    
    ina = P*ghk(v, nai, nao,1)*Area * g_factor	
    ik = P*ghk(v, ki, ko,1)*Area * g_factor
    iampa = ik + ina
    
}

DERIVATIVE state {
    lr=eta(cai)
    
    A' = -A/taur
    B' = -B/taud
    
    U_SE = U_SE_init * (1 + U_SE_factor)
    x' = z / tau_rec - U_SE * x * I
    y_rel' = -y_rel / tau_in + U_SE * x * I
    z = 1 - x - y_rel
    gluti' = -gluti / tau_takeover + y_rel * glut_factor 
    
}

FUNCTION ghk(v(mV), ci(mM), co(mM),z) (0.001 coul/cm3) {
    LOCAL arg, eci, eco
    arg = (0.001)*z*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(arg)
    eci = ci*efun(-arg)
    ghk = (0.001)*z*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
	efun = 1 - z/2
    }else{
	efun = z/(exp(z) - 1)
    }
}

FUNCTION eta(ci (mM)) { 
    LOCAL inv, P1, P2, P3, P4
    P1=100	
    P2=P1*1e-4	
    P4=1e3
    P3=3		
    
    ci=(ci-1e-4)*1e3 	
    
    inv=P4 + P1/(P2+ci*ci*ci) 
    eta=1/inv
}	

FUNCTION Omega(ci (mM)) {
    ci=(ci-1e-4)*1e3	
    Omega=0.25+1/(1+exp(-(ci-alpha2)*beta2))-0.25/(1+exp(-(ci-alpha1)*beta1))
}

NET_RECEIVE(weight (uS)) { 	
    
    if( flag==0 ) {  
	I = 10
	net_send(0.1,2)
	}	
    if( flag==2 ) { 
	I = 0
	
    }

    
    
    
    
    A = A + Pmax * factor * gluti
    B = B + Pmax * factor * gluti
    
}