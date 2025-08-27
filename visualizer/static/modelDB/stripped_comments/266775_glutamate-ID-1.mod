NEURON {
	POINT_PROCESS glutamate
	RANGE tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda
	RANGE erev, g, i
	RANGE i_ampa, i_nmda, g_ampa, g_nmda, ratio, I, G, mg, q, block, alpha, beta
	RANGE ampa_scale_factor, nmda_scale_factor
    RANGE damod, maxModNMDA,max2NMDA,maxModAMPA,max2AMPA,l1NMDA,l2NMDA,l1AMPA,l2AMPA
	
	NONSPECIFIC_CURRENT i
	USEION cal WRITE ical VALENCE 2
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	erev        = 0.0       (mV)
	
	tau1_ampa   = 1.9       (ms)
    tau2_ampa   = 4.8       (ms)  
    tau1_nmda   = 5.52      (ms)  
    tau2_nmda   = 231       (ms)  
    
    ratio       = 1         (1)   
    mg          = 1         (mM)
    alpha       = 0.062
    beta        = 3.57
    q           = 2               
    
    nmda_scale_factor = 1
    ampa_scale_factor = 1
    
    ca_ratio_ampa = 0.005
    ca_ratio_nmda = 0.1
    
    maxModNMDA  = 1
    max2NMDA    = 1
    maxModAMPA  = 1
    max2AMPA    = 1
    damod       = 0
    l1NMDA      = 0
    l2NMDA      = 0
    l1AMPA      = 0
    l2AMPA      = 0
}


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor_nmda
	factor_ampa
	i_ampa
	i_nmda
	g_ampa
	g_nmda
	block
	I
	G
	ical (nA)
}


STATE {
	A (uS)
	B (uS)
	C (uS)
	D (uS)
}



INITIAL {
	LOCAL tp
	if (tau1_nmda/tau2_nmda > .9999) {
		tau1_nmda = .9999*tau2_nmda
	}
	if (tau1_ampa/tau2_ampa > .9999) {
		tau1_ampa = .9999*tau2_ampa
	}
	
	
	A           = 0
	B           = 0
	tp          = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)
	factor_nmda = -exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda)
	factor_nmda = 1/factor_nmda
	
	
	C           = 0
	D           = 0
	tp          = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)
	factor_ampa = -exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa)
	factor_ampa = 1/factor_ampa
}




BREAKPOINT {
	SOLVE state METHOD cnexp
	
	
	g_nmda = (B - A) * modulation(maxModNMDA,max2NMDA,l1NMDA,l2NMDA)
	block  = MgBlock()
	i_nmda = g_nmda * (v - erev) * block * nmda_scale_factor
	
	
	g_ampa = (D - C) * modulation(maxModAMPA,max2AMPA,l1AMPA,l2AMPA)
	i_ampa = g_ampa * (v - erev) * ampa_scale_factor
	
	
	G = g_ampa + g_nmda
	I = i_ampa + i_nmda
	
	
	ical = i_ampa*ca_ratio_ampa  + i_nmda*ca_ratio_nmda
    i = i_ampa*(1-ca_ratio_ampa) + i_nmda*(1-ca_ratio_nmda)
}



DERIVATIVE state {
	A' = -A/tau1_nmda*q
	B' = -B/tau2_nmda*q
	C' = -C/tau1_ampa*q
	D' = -D/tau2_ampa*q
}



NET_RECEIVE(weight (uS)) {
	A = A + weight*factor_nmda
	B = B + weight*factor_nmda
	C = C + weight*factor_ampa*ratio
	D = D + weight*factor_ampa*ratio
}



FUNCTION MgBlock() {
    
    MgBlock = 1 / (1 + mg * exp(-alpha * v) / beta )
    
}

FUNCTION modulation(m1,m2,l1,l2) {
    
    
    modulation = 1 + damod * ( (m1-1)*l1 + (m2-1)*l2 )
    if (modulation < 0) {
        modulation = 0
    } 
}