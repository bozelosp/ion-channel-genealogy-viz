UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (mW) = (milliwatt)
    (photons)  = (1)
}

NEURON { 
	SUFFIX chanrhod
	
    NONSPECIFIC_CURRENT	icat
    
    RANGE photons, flux, irradiance, U, U0, U1
    
    RANGE channel_density, gdens1, gdens2
    
    RANGE x, y, z
    
	
	
    GLOBAL source_photons, source_flux, source_irradiance
    
    RANGE Ka1, Ka2
    GLOBAL Kd1, Kd2
    GLOBAL Kr
    RANGE e12, e21, e12dark, e21dark
    RANGE delta1, delta2
    RANGE o10, o20
    RANGE Islow, Ifast
    RANGE c_1, c_2, b
    GLOBAL h, c
    RANGE amp, photon_energy
    RANGE wavelength
    RANGE phi0, phio, gcat    
    GLOBAL sigma_retinal, gcat1, gcat2, ecat, Imax, gamma, tChR
    
    GLOBAL epsilon1, epsilon2
    
    RANGE Tx, phi
    
    GLOBAL tstimon, tstimoff
}



PARAMETER {
    
    channel_density        = 1.3e10              (1/cm2) 
    
    

    
    
    gcat1    = 50e-15    (mho)   
    gcat2    = 250e-17   (mho)   
    
    
    
    

    gamma=0.05 

    sigma_retinal = 1.2e-16  (cm2)        

    epsilon1 = 0.5             (1)         
	
	epsilon2 = 0.1                      
    ecat     = 0      (mV)     
                               
                               

    Tx      = 1       (1)      
    vshift  = 0        (mV)     

	
	tChR=1.3 (ms) 

    x = 0 (1) 
	y = 0 (1)
	z = 0 (1)
	Ka1 = 0.5 
	Ka2 = 0.12 

	
	
	Kd1 = 0.13 
	Kd2 = 0.025 

	Kr = 0.0004 

	e12 = 0.053 
	e21 = 0.023 
	e12dark = 0.022
	e21dark = 0.011

	delta1= 0.03 (ms)
	delta2= 0.15 (ms)

	h = 6.6260693e-34           (m2 kg/s)  
    c = 299792458.0             (m/s)      
	wavelength = 4.45e-7
}

ASSIGNED {  
    v           (mV)
    icat        (mA/cm2)
    gdens1        (mho/cm2)
    gdens2        (mho/cm2)
    source_irradiance  (W/cm2)      
    source_photons     (photons/ms) 
    source_flux        (photons/ms cm2) 
    irradiance          (W/cm2)      
    flux               (photons/ms cm2) 
    phi                (photons/ms) 
    U
    U0
    U1
    Imax
    Islow
    Ifast
    c_1
    c_2
    b
    amp
    photon_energy
    phi0
    phio
    gcat
    tstimon
    tstimoff
}

STATE { 
	o1 o2 c1 c2
}

INITIAL {
    irradiance = 0
    flux = 0
    phi = 0
    Islow=0
    Ifast=0
    
	
    tstimon = 0
    tstimoff = 0

    
	c1 = 1 
	c2 = 0
	o1 = 0
	o2 = 0

	phio = 0
    o10=0
    o20=0
}

BREAKPOINT {
    irradiance = source_irradiance * Tx 
    flux      = source_flux * Tx           
    phi       = flux             * sigma_retinal
                
                

	U=v-ecat-vshift-75
	U0=40
	U1=15
	Imax=(v-ecat-vshift)*gcat1*channel_density 

    b= (Kd1+Kd2+e12dark+e21dark)/2

    c_1 = 0.1029797709
    c_2 = 0.0398631371

    gcat=(o1+gamma*o2)*gcat1*channel_density

	if (phi>0) {
		Ka1 = epsilon1 * phi * (1 - exp( -(t - tstimon) / tChR)) 
		Ka2 = epsilon2 * phi * (1 - exp( -(t - tstimon) / tChR))
		
		
		e12=0.053
		e21=0.023
		
   		
		icat = Imax * (o1 + gamma * o2) * (1-exp(-U/U0))/(U/U1) 
		o10=o1
		o20=o2
		phio=phi
	} else {
		Ka1 = epsilon1 * phio * (exp (- ((t-tstimoff) / tChR)) - exp( -(t - tstimon) / tChR) ) 
		Ka2 = epsilon2 * phio * (exp (- ((t-tstimoff) / tChR)) - exp( -(t - tstimon) / tChR) ) 
		e12 = 0.022
		e21 = 0.011
		
    	Islow = Imax * ( ( (delta2 - (Kd1 + (1 - gamma) * e12dark)) * o10 + (( 1 - gamma) * e21dark + gamma*(delta2-Kd2)) * o20 ) / (delta2 - delta1) )
    	Ifast = Imax * ( ( (Kd1 + (1 - gamma) * e12dark - delta1) * o10 + (-(1 - gamma)*e21dark+gamma*(Kd2-delta1))*o20 ) / (delta2 - delta1) )
    	icat = Islow * exp(-delta1*(t-tstimoff)) + Ifast*exp(-delta2*(t-tstimoff)) 
	}

    SOLVE states METHOD cnexp
    if (o1>1){o1=1}
    if (o1<0){o1=0}
	if (o2>1){o2=1}
    if (o2<0){o2=0}
    if (c1>1){c1=1}
    if (c1<0){c1=0}
    if (c2>1){c2=1}
    if (c2<0){c2=0}
    c1    = 1 - o1 - o2 - c2
}

DERIVATIVE states {  
	o1' = Ka1*c1 - (Kd1 + e12)*o1 + e21*o2
	o2' = Ka2*c2 + e12*o1 - (Kd2+e21)*o2
	
	c2' = Kd2*o2 - (Ka2+Kr)*c2
}