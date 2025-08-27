NEURON { 
	SUFFIX nadifus
	
	USEION na READ nao, nai, ina WRITE nai
	RANGE inabar
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
	nao		(mM)
	ina		(mA/cm2)
	k1buf		(/mM-ms)
	k2buf		(/ms)
	inabar		(mA/cm2)
	nai0 = 12	(mM)		
}

ASSIGNED {
	nai		(mM)
	vol[NANN]	(1)	
}

STATE {
	na[NANN]	(mM) 
	NaBuffer[NANN]  (mM)
	Buffer[NANN]    (mM)
}

BREAKPOINT {
	SOLVE state METHOD sparse
	ina = inabar
}

LOCAL coord_done

INITIAL {
	if (coord_done == 0) {
		coord_done = 1
		coord()
	}
	
	
	FROM i=0 TO NANN-1 {
		na[i] = nai0
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
	COMPARTMENT i, diam*diam*vol[i] {na NaBuffer Buffer}
	LONGITUDINAL_DIFFUSION j, DFree*diam*diam*vol[j] {na}
	~ na[0] << (-ina*PI*diam*frat[0]/(FARADAY))
	FROM i=0 TO NANN-2 {
		~ na[i] <-> na[i+1] (DFree*frat[i+1], DFree*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ na[i] + Buffer[i] <-> NaBuffer[i] (k1buf*dsqvol,k2buf*dsqvol)
	}
	
	nai = (na[0] + na[1] + na[2] + na[3])/4
}