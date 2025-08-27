NEURON {
	SUFFIX kcl
	USEION k  READ ki , ko  WRITE ik, ek   VALENCE 1
	USEION cl READ cli, clo WRITE icl, ecl VALENCE -1
	RANGE veq, gsum
	RANGE ik, icl, Pratio, Pk ,Pcl
	
	USEION na READ nai, nao WRITE ena      VALENCE 1
}

UNITS {
	(mV)     = (millivolt)
	(mA)     = (milliamp)
	(S)      = (siemens)
 	(molar)  = (1/liter)	
	(mM)     = (millimolar)
    FARADAY  = 9.6485e4 (coulombs)
 	R        = 8.3145   (volt coul/kelvin) 	
 	(um)     = (micrometer)
}

PARAMETER {
	gsum (S/cm2)          
	veq (mV)
		
	zk  = 1  (1)         
	zcl = -1 (1)
	zna = 1  (1)	

	rnak = 1 
	
	a = 150
	b = 5
	c = 3
}

ASSIGNED {
    ki  (mM)
    ko  (mM)
    cli (mM)
    clo (mM)
    nai (mM)
    nao (mM)
    
    ik    (mA/cm2)    
    icl   (mA/cm2) 
	    
    Pratio (1)
    
    ek  (mV)
    ecl (mV)
    ena (mV)
	
    Pcl (cm/sec)
    Pk  (cm/sec)
    
    Gcl (sec S / cm3)
    Gk  (sec S / cm3)
    
    v       (mV) 
    celsius (degC)
    diam    (um)
	      
    T     (kelvin)
    E     (volt)

    beta  (1)
    vt (mV) 
}



INITIAL {
    T = kelvinfkt(celsius)
    vt = R*T/FARADAY      
    ek  = NERNST(ko, ki  , zk)   
    ecl = NERNST(clo, cli, zcl)
    ena = NERNST(nao, nai, zna)
    
    
    beta = exp(veq*(1e-3) * FARADAY / (R*T) ) 
    Pratio =  ( cli - beta*clo ) / ( beta*ki - ko )
    Gcl =          GOLDMAN(veq*(1e-3), cli,clo,zcl) / (veq-ecl)  
    Gk  = Pratio * GOLDMAN(veq*(1e-3), ki,ko,zk)    / (veq-ek) 
    Pcl = gsum / (Gcl+Gk)
    Pk  = Pratio * Pcl
    
    if (Pratio < 0) {
        
        Pratio = 0
        Pk = 0
        Pcl = (1e9) * gsum / (  (FARADAY*FARADAY*FARADAY*veq / (R*R*T*T) ) * ( (clo*cli)/(cli-clo) )  )
    }
}

BREAKPOINT {
    T = kelvinfkt(celsius)
    vt = R*T/FARADAY      
    E = v*(1e-3)
    ek  = NERNST(ko, ki  , zk)   
    ecl = NERNST(clo, cli, zcl)
    ena = NERNST(nao, nai, zna)
    
    
    beta = exp(veq*(1e-3) * FARADAY / (R*T) )
    Pratio = ( cli - beta*clo ) / ( beta*ki - ko )
    Gcl =          GOLDMAN(veq*(1e-3),cli,clo,zcl) / (veq-ecl) 
    Gk  = Pratio * GOLDMAN(veq*(1e-3),ki,ko,zk)    / (veq-ek) 
    Pcl = gsum / (Gk+Gcl)
    Pk  = Pratio * Pcl
    
    if (Pratio < 0) {
        
        Pratio = 0
        Pk  = 0
        Pcl =  (1e9) *gsum / (  (FARADAY*FARADAY*FARADAY*veq / (R*R*T*T) ) * ( (clo*cli)/(cli-clo) )  ) 
    }
    
    ik  = Pk  * GOLDMAN(E, ki, ko, zk)
    icl = Pcl * GOLDMAN(E, cli, clo, zcl)
    
}

FUNCTION GOLDMAN (V (volt), ci (mM), co (mM), z (1)) (mA sec/cm3) {
    LOCAL zeta, ezeta
 	zeta = z*V/vt
    ezeta = exp(-zeta)
 	GOLDMAN = (z*FARADAY*zeta) * (1e-3) * (ci - co * ezeta )/(1-ezeta) 
}

FUNCTION kelvinfkt (t (degC) ) (kelvin) {
 	kelvinfkt = 273.19 + t
} 

FUNCTION NERNST (co (mM), ci (mM), z (1)) (mV) {
    NERNST = (1e3) * vt/z * log(co/ci)
}