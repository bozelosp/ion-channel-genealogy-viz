NEURON {
	 SUFFIX ip3dif
 	 USEION ip3 READ ip3i WRITE ip3i VALENCE 1
  	 GLOBAL vol, DIP3, ip3i0, kdegr
     RANGE ip3i
	 THREADSAFE
}

DEFINE NANN 12

UNITS {
  	(molar) = (1/liter)
  	(mM)    = (millimolar)
  	(uM)    = (micromolar)
  	(um)    = (micron)
  	(mA)    = (milliamp)
  	FARADAY = (faraday)  (coulomb)
  	PI      = (pi)       (1)
}

PARAMETER {
  	kdegr = 0.14e-3 (/ms)  
  	DIP3 = 0.283(um2/ms)
  	ip3i0 = 0.16e-3 (mM)   
}


ASSIGNED {
  	diam      (um)
 	ip3i      (mM)
  	vol[NANN]  		
					
					
}

STATE {
  	ip3[NANN]       (mM) <1e-6>
}

LOCAL factors_done

BREAKPOINT { SOLVE state METHOD sparse}


INITIAL {
   	if (factors_done == 0) {   
     	factors_done = 1       
      	factors()              
       }

  	ip3i = ip3i0
  	FROM i=0 TO NANN-1 {
    	ip3[i] = ip3i
	}
}


LOCAL frat[NANN]  

PROCEDURE factors() {
  	LOCAL r, dr2
  	r = 1/2                
  	dr2 = r/(NANN-1)/2     
						   
  	vol[0] = 0
  	frat[0] = 2*r
 	 FROM i=0 TO NANN-2 {
    		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2  
    		r = r - dr2
    		frat[i+1] = 2*PI*r/(2*dr2)  
										
   		 r = r - dr2
    		vol[i+1] = PI*(r+dr2/2)*2*dr2  
  	}
}

LOCAL dsq, dsqvol  

KINETIC state {
  	COMPARTMENT i, diam*diam*vol[i]*0.81 {ip3 ip3i0}  
	
	dsq = diam*diam
	
  	FROM i=0 TO NANN-2 {
   		 ~ ip3[i] <-> ip3[i+1]  (DIP3*frat[i+1], DIP3*frat[i+1])
  	}

  	
  	FROM i=0 TO NANN-1 {
    		dsqvol = dsq*vol[i]*0.81 						 
    		~ ip3[i] <-> ip3i0  (kdegr*dsqvol, kdegr*dsqvol) 
  	}

  	ip3i = ip3[0]
}