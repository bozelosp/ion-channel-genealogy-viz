NEURON {
	  SUFFIX caextscale						
	  USEION ca READ ica WRITE cao		
	  GLOBAL cabath						
	  RANGE fhspace, txfer, Vol_peri, L, nseg, SA, Volratio, lseg		
	}


	UNITS {
	  (mV)    = (millivolt)
	  (mA)    = (milliamp)
	  FARADAY = 96500 (coulombs)
	  (molar) = (1/liter)
	  (mM)    = (millimolar)
	  (um)	  = (micrometer)
	  PI      = (pi) (1)
	}


	PARAMETER {
	  cabath   =  2 (mM)        	
	  fhspace = .03 (um)  			
	  txfer   =  4511.0 (ms)  		
	  
	  VolSchild = 1.41372E-08 (cm3)
	  Vol_periSchild = 1.46136E-09 (cm3)	
	  
	}


	ASSIGNED { 
	
		ica  		(mA/cm2)
		SA 			(cm2)
		Vol_peri	(cm3)
		diam		(um)
		nseg		(1)
		L			(um)
		lseg		(cm)
		Vol			(cm3)
		Volratio	(1)

	}


	STATE { cao  (mM) }


	BREAKPOINT { SOLVE state METHOD derivimplicit }


	INITIAL {
		lseg=(1e-4)*L/nseg
		SA = PI*(1e-4)*diam*lseg
		Vol = (PI*((1e-4)*(diam/2))^2*lseg)
		Volratio = Vol_periSchild/VolSchild
		Vol_peri = (PI*((1e-4)*((diam+fhspace)/2))^2*lseg)-Vol
		
	}


	DERIVATIVE state {
	  cao' = ica*SA/(2*Vol_peri*FARADAY) + (cabath - cao)/txfer 
	}