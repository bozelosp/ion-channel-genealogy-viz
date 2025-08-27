NEURON {
	SUFFIX cadifus
	USEION ca READ cao, cai, ica WRITE cai, ica
	RANGE icabar
	GLOBAL vol, Buffer0
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
	cao		(mM)
	ica		(mA/cm2)
	k1buf		(/mM-ms)
	k2buf		(/ms)
	icabar		(mA/cm2)
}

ASSIGNED {
	cai		(mM)
	vol[NANN]	(1)	
}

STATE {
	ca[NANN]	(mM) 
	CaBuffer[NANN]  (mM)
	Buffer[NANN]    (mM)
}

BREAKPOINT {
	SOLVE state METHOD sparse
	ica = icabar
}

LOCAL coord_done

INITIAL {
	if (coord_done == 0) {
		coord_done = 1
		coord()
	}
	
	
	FROM i=0 TO NANN-1 {
		ca[i] = cai
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
	COMPARTMENT i, diam*diam*vol[i] {ca CaBuffer Buffer}
	LONGITUDINAL_DIFFUSION j, DFree*diam*diam*vol[j] {ca}
	~ ca[0] << (-ica*PI*diam*frat[0]/(2*FARADAY))
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DFree*frat[i+1], DFree*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol,k2buf*dsqvol)
	}
	cai = ca[0]
}