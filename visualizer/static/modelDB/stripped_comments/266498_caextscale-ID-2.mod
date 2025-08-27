NEURON {
	  SUFFIX caextscale						
	  USEION ca READ ica WRITE cao		
	  GLOBAL cabath						
	  RANGE fhspace, txfer, Vol_peri, L, nseg, SA, lseg		
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
	  fhspace = 1 (um)  			
	  txfer   =  4511.0 (ms)  		
  
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

	}


	STATE { cao  (mM) }


	BREAKPOINT { SOLVE state METHOD derivimplicit }


	INITIAL {
		lseg=(1e-4)*L/nseg
		SA = PI*(1e-4)*diam*lseg
		Vol = (PI*((1e-4)*(diam/2))^2*lseg)
		Vol_peri = (PI*((1e-4)*((diam+fhspace)/2))^2*lseg)-Vol
	}


	DERIVATIVE state {
	  cao' = ica*SA/(2*Vol_peri*FARADAY) + (cabath - cao)/txfer 
	}