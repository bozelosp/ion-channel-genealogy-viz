NEURON {
	SUFFIX hcn
	NONSPECIFIC_CURRENT i
	RANGE i, ehcn, g, gbar
    EXTERNAL apc_metap, fpc_metap
	GLOBAL a0, b0, ah, bh, ac, bc, aa0, ba0
	GLOBAL aa0, ba0, aah, bah, aac, bac
	GLOBAL kon, koff, b, bf, gca, shift
    GLOBAL Vhalf, vh1, vh2, vh_shift, avh1, avh2, avh_shift
	RANGE ai
}

UNITS {
	(mV)	= (millivolt)
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA)	= (milliamp)
	(S)	= (siemens)
}

PARAMETER {
	gbar    = 0.0		(S/cm2)
    ehcn    = -44       (mV) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    a0      = .00032743		(/ms)		
    b0      = .00029334		(/ms)
    ac      = -0.1103		(/mV)
    bc      = 0.1025	    	(/mV)

    aa0     = 0.0011    	(/ms)		
    ba0     = 0.0164 		(/ms)
    aac     = -0.0774		(/mV)
    bac     = 0.1486   		(/mV)

    
    Vhalf = -90.0           (mV)
    vh1   = 1.057
    vh2   = 79.3
    avh1   = 1.886        
    avh2 = 164.1       

    
    kon     = 3.086		    (/mM-ms)	
    koff    = 4.4857e-05	(/ms)
    b       = 80
    bf      = 8.94

	ai	= 1e-05		(mM)		
	gca     = 1				
	shift   = 0		(mV)		
	q10v    = 4				
	q10a    = 1.5				
	celsius			(degC)
}

ASSIGNED {
	v	(mV)
	g	(S/cm2)
	i	(mA/cm2)
	alpha	(/ms)
	beta    (/ms)
	alphaa	(/ms)
	betaa	(/ms)

    vh_shift
    avh_shift

    ah  (mV)
    bh  (mV)
    aah (mV)
    bah (mV)
}

STATE {
	c
	cac
	o
	cao
}

INITIAL {
    setVhalf(Vhalf)
    SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*(o + cao*gca)
	i = g*(v-ehcn)
}

KINETIC kin {
	LOCAL qa
	qa = q10a^((celsius-22 (degC))/10 (degC)) 
	
	rates(v)
	~ c <-> o       (alpha, beta)
	~ c <-> cac     (kon*qa*ai/bf,koff*qa*b/bf)
	~ o <-> cao     (kon*qa*ai, koff*qa)
	~ cac <-> cao   (alphaa, betaa)
	CONSERVE c + cac + o + cao = 1
}

PROCEDURE rates(v(mV)) {
	LOCAL qv
	qv = q10v^((celsius-22 (degC))/10 (degC)) 
	
	if (v > -200) {
		alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))
	} else {
		alpha = a0*qv / (1 + exp(-((-200)-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-((-200)-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-((-200)-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-((-200)-bah-shift)*bac))
	}
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vh_shift = Vhalf*vh1+vh2
    
    avh_shift = Vhalf*avh1+avh2

    ah      = -87.7 + vh_shift
    bh      = -51.7 + vh_shift
    aah     = -94.2 + avh_shift
    bah     = -35.5 + avh_shift
}