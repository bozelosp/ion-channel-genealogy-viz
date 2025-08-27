NEURON {
	SUFFIX kdifus
	
	USEION k READ ko, ki, ik WRITE ki
	RANGE ikbar
	GLOBAL vol 
}

DEFINE NANN  4

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	FARADAY = (faraday)	 (10000 coulomb)
	PI	= (pi) (1)
}

PARAMETER {
	DFree = .6	(um2/ms)
	diam		(um)
	ko		(mM)
	ik		(mA/cm2)
	k1buf		(/mM-ms)
	k2buf		(/ms)
	ikbar		(mA/cm2)
	ki0 = 150 	(mM)		
}

ASSIGNED {
	ki		(mM)
	vol[NANN]	(1)	
}

STATE {
	k[NANN]	(mM) 	
	KBuffer[NANN]  (mM)
	Buffer[NANN]    (mM)
}

BREAKPOINT {
	SOLVE state METHOD sparse
	ik = ikbar
}

LOCAL coord_done

INITIAL {
	if (coord_done == 0) {
		coord_done = 1
		coord()
	}
	
	
	FROM i=0 TO NANN-1 {
		k[i] = ki0
	}
}

LOCAL frat[NANN] 	

PROCEDURE coord() {
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
	COMPARTMENT i, diam*diam*vol[i] {k KBuffer Buffer}
	LONGITUDINAL_DIFFUSION j, DFree*diam*diam*vol[j] {k}
	~ k[0] << (-ik*PI*diam*frat[0]/(FARADAY))
	FROM i=0 TO NANN-2 {
		~ k[i] <-> k[i+1] (DFree*frat[i+1], DFree*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ k[i] + Buffer[i] <-> KBuffer[i] (k1buf*dsqvol,k2buf*dsqvol)
	}
	
	ki = (k[0] + k[1] + k[2] + k[3])/4
}