NEURON {
    POINT_PROCESS ghknmda
    USEION na WRITE ina
    USEION k WRITE ik
    USEION ca READ cai, cao WRITE ica
    USEION glut READ gluti VALENCE 0
    
    RANGE taur, taud
    RANGE inmda
    
    RANGE P, mg, Pmax
    RANGE  mgb, ica, Area, mgb_k, mg_ref
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
    taur=5 (ms) <1e-9,1e9>
    taud = 50 (ms) <1e-9,1e9>
    cai = 100e-6(mM)	
    cao = 2		(mM)
    nai = 18	(mM)	
    nao = 140	(mM)
    ki = 140	(mM)	
    ko = 5		(mM)
    celsius		(degC)
    mg = 1		(mM) 
    Pmax=1e-6   (cm/s)	
    
    Area = 1 (cm2)
    mgb_k = 0.062 (/mV)
    mg_ref = 3.57 (mM)
}

ASSIGNED {
    ina     (nA)
    ik      (nA)
    ica     (nA)
    v (mV)
    P (cm/s)
    factor
    mgb	(1)
    inmda	(nA)
    
    
    gluti (mM)
}

STATE {
    A (cm/s)
    B (cm/s)
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
    
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    P=B-A
    mgb = mgblock(v)
    
    
    
    ina = P*mgb*ghk(v, nai, nao,1)*Area	
    ica = P*10.6*mgb*ghk(v, cai, cao,2)*Area
    ik = P*mgb*ghk(v, ki, ko,1)*Area
    inmda = ica + ik + ina
    
    
}

DERIVATIVE state {
    A' = -A/taur
    B' = -B/taud
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

FUNCTION mgblock(v(mV)) (1){
    TABLE 
    DEPEND mg
    FROM -140 TO 80 WITH 1000 
    
    
    
    mgblock = 1 / (1 + exp(mgb_k * -v) * (mg / mg_ref)) 
    
    
    
    
}

NET_RECEIVE(weight (uS)) { 	
    
    
    
    if (flag == 0 ) {
	net_send(.01,2)
    }
    if (flag == 2 ) {
	A = A + Pmax*factor * gluti
	B = B + Pmax*factor * gluti
    }
}