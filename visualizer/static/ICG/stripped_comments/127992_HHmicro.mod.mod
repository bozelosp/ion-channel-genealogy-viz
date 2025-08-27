INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micron)
} 


NEURON {
    SUFFIX HHmicro
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE el
    RANGE gnabar, gkbar, gna, gk, gl
    RANGE m0h0,m1h0,m2h0,m3h0,m0h1,m1h1,m2h1,m3h1
    RANGE n0, n1, n2, n3, n4
    RANGE gamma_na, gamma_k
    RANGE Nna, Nk
    GLOBAL m_inf, h_inf, n_inf, tau_m, tau_h, tau_n
    RANGE seed
    THREADSAFE 
} 


PARAMETER {
    gnabar = 0.0   (S/cm2)     
    gkbar  = 0.036  (S/cm2)     
    gl     = 0.0 (S/cm2)     
    el       = -54.3 (mV)       
    
    gamma_na = 10  (pS)		
    gamma_k  = 10  (pS)		
    seed = 5061983              
} 


STATE {
    m h n
} 


ASSIGNED {
    ina   (mA/cm2)
    ik    (mA/cm2)
    il    (mA/cm2)
    gna   (S/cm2)
    gk    (S/cm2)
    ena  (mV)
    ek   (mV)
    
    dt    (ms)
    area  (um2)
    celsius  (degC)
    v  (mV)
    
    Nna			 	
    Nk			 	

    m0h0			
    m1h0			
    m2h0			
    m3h0			
    m0h1			
    m1h1			
    m2h1			
    m3h1			

    n1				
    n2				
    n3				
    n4				
    n5				
    
    m_inf h_inf n_inf
    tau_m (ms) tau_h (ms) tau_n (ms)
    
} 



INITIAL {
    m = 0
    h = 0
    n = 0
    
    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   
    
    VERBATIM
		
    ENDVERBATIM

    m0h0 = 0
    m1h0 = 0
    m2h0 = 0
    m3h0 = 0
    m0h1 = Nna		
    m1h1 = 0
    m2h1 = 0
    m3h1 = 0

    n1   = Nk		
    n2   = 0
    n3   = 0
    n4   = 0
    n5   = 0
    
    set_seed(seed)
    
} 


BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = m3h1*((1e-12)*gamma_na)/((1e-8)*area)
    gk = n5*((1e-12)*gamma_k)/((1e-8)*area)
    ina = gna * (v-ena)
    ik = gk * (v-ek)
    il = gl * (v-el)
} 


DERIVATIVE states {   
    UNITSOFF
    m' = m
    h' = h
    n' = n
    UNITSON
    noise()
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


PROCEDURE noise() {
    LOCAL p,am,bm,ah,bh,an,bn,m0h0_merk,m1h0_merk,m2h0_merk,m3h0_merk,m0h1_merk,m1h1_merk,m2h1_merk,m3h1_merk,n1_merk,n2_merk,n3_merk,n4_merk,n5_merk,rnd,den
    
    
    am = alpham(v)
    bm = betam(v)
    
    
    ah = alphah(v)
    bh = betah(v)
    
    
    an = alphan(v)
    bn = betan(v)

    m0h0_merk = m0h0
    m1h0_merk = m1h0
    m2h0_merk = m2h0
    m3h0_merk = m3h0
    m0h1_merk = m0h1
    m1h1_merk = m1h1
    m2h1_merk = m2h1
    m3h1_merk = m3h1

    n1_merk = n1
    n2_merk = n2
    n3_merk = n3
    n4_merk = n4
    n5_merk = n5
    
    
    
    
    den = ah+3*am
    p = 1-exp(-dt*den)
    FROM ii=1 TO m0h0_merk {
	
	if (scop_random() <= p)	{					
	    if (scop_random() <= ah/den) {       		        
		m0h0=m0h0-1
		m0h1=m0h1+1
	    }
	    else {							
		m0h0=m0h0-1
		m1h0=m1h0+1
	    }
	}
    }
    
    
    den = ah+2*am+bm
    p = 1-exp(-dt*den)
    FROM ii=1 TO m1h0_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= ah/den) {	
		m1h0=m1h0-1
		m1h1=m1h1+1
	    }
	    else if (rnd <= (ah+2*am)/den) {
		m1h0=m1h0-1
		m2h0=m2h0+1
	    }
	    else {
		m1h0=m1h0-1
		m0h0=m0h0+1
	    }
	}
    }
    
    
    den = ah+am+2*bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m2h0_merk {
	if (scop_random() <= p){
	    rnd = scop_random()
	    if (rnd <= ah/den)	{
		m2h0=m2h0-1
		m2h1=m2h1+1
	    }
	    else if (rnd <= (ah+am)/den) {
		m2h0=m2h0-1
		m3h0=m3h0+1
	    }
	    else {
		m2h0=m2h0-1
		m1h0=m1h0+1
	    }
	}
    }
    
    
    den = ah+3*bm
    p = 1-exp(-dt*den)
    FROM ii=1 TO m3h0_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= ah/den) {	
		m3h0=m3h0-1
		m3h1=m3h1+1
	    }
	    else {	
		m3h0=m3h0-1
		m2h0=m2h0+1
	    }
	}
    }
    
    
    
    
    
    den = bh+3*am
    p=1-exp(-dt*den) 
    FROM ii=1 TO m0h1_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m0h1=m0h1-1
		m0h0=m0h0+1
	    }
	    else {
		m0h1=m0h1-1
		m1h1=m1h1+1
	    }
	}
    }
    
    
    den = bh+2*am+bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m1h1_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m1h1=m1h1-1
		m1h0=m1h0+1
	    }
	    else if (rnd <= (bh+2*am)/den) {	
		m1h1=m1h1-1
		m2h1=m2h1+1
	    }
	    else {
		m1h1=m1h1-1
		m0h1=m0h1+1
	    }
	}
    }
    
    
    den = bh+am+2*bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m2h1_merk {
	if (scop_random()<= p){
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m2h1=m2h1-1
		m2h0=m2h0+1
	    }
	    else if (rnd <= (bh+am)/den) {	
		m2h1=m2h1-1
		m3h1=m3h1+1
	    }
	    else {
		m2h1=m2h1-1
		m1h1=m1h1+1
	    }
	}
    }
    
    
    den = bh+3*bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m3h1_merk {
	if (scop_random()<= p) {
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m3h1=m3h1-1
		m3h0=m3h0+1
	    }	
	    else {
		m3h1=m3h1-1
		m2h1=m2h1+1
	    }
	}
    }
    
    
    
    
    
    den = 4*an
    p = 1-exp(-dt*den) 
    FROM ii=1 TO n1_merk {
	if (scop_random() <= p) {
	    n1=n1-1
	    n2=n2+1	
	}
    }
    
    
    
    den = 3*an+bn
    p = 1-exp(-dt*den) 
    FROM ii=1 TO n2_merk {
	if (scop_random() <= p) {
	    if (scop_random() <= 3*an/den) {	
		n2=n2-1
		n3=n3+1
	    }
	    else {
		n2=n2-1
		n1=n1+1
	    }
	}
    }
    
    
    
    den = 2*an+2*bn
    p = 1-exp(-dt*den)
    FROM ii=1 TO n3_merk {
	if (scop_random() <= p) {
	    if (scop_random() <= 2*an/den) {
		n3=n3-1
		n4=n4+1
	    }
	    else {
		n3=n3-1
		n2=n2+1
	    }
	}
    }
    
    
    
    den = an+3*bn
    p = 1-exp(-dt*den)
    FROM ii=1 TO n4_merk {
	if (scop_random() <= p) {
	    if (scop_random() <= an/den) {
		n4=n4-1
		n5=n5+1
	    }
	    else {
		n4=n4-1
		n3=n3+1
	    }
	}
    }
    
    
    
    den = 4*bn
    p = 1-exp(-dt*den)
    FROM ii=1 TO n5_merk {
	if (scop_random() <= p) {
	    n5=n5-1
	    n4=n4+1	
	}
    }
}