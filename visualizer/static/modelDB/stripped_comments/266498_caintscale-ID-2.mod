NEURON {
	  SUFFIX caintscale					
	  USEION ca READ ica WRITE cai		
	  GLOBAL nb							
	  RANGE a, ku, kr, Bi, diffOc, SA, Vol, L, nseg	
	}


	UNITS {
	  (mV)   	= (millivolt)
	  (um)		= (micrometer)
	  (mA)   	= (milliamp)
	  FARADAY 	= 96500 (coulombs)
	  (molar)	= (1/liter)
	  (mM)    	= (millimolar)
	  PI 		= (pi) (1)
	}


	PARAMETER {
	  a = 15e-4 (cm)  			
	  ku = 100 (/mM/ms) 		
	  kr = 0.238 (/ms)		 	
	  nb = 4 					
	  Bi = 0.001 (mM) 			
	  
	  
	  
	  
	}


	ASSIGNED { 
		
		
		ica  	(mA/cm2)
		L 		(um)
		nseg	(1)
		lseg	(cm)
		SA		(cm2)
		Vol		(cm3)
		diam	(um)
		
		diffOc (/ms)


	}


	STATE { cai  (mM) Oc } 


	BREAKPOINT { SOLVE state METHOD derivimplicit }


	INITIAL {

	Oc=.05
	lseg=(1e-4)*L/nseg
	SA = PI*(1e-4)*diam*lseg
	Vol = (PI*((1e-4)*(diam/2))^2*lseg)

	}


	DERIVATIVE state {
	  LOCAL diffOc											
	  diffOc=ku*cai*(1-Oc)-kr*Oc 							
	  Oc'=diffOc											
	  cai' = -ica*(SA)/Vol/(2*FARADAY) - (nb*Bi*diffOc)		
	}