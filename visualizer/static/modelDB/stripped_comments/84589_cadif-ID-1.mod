NEURON {
	SUFFIX cadifus
	USEION ca READ cai, ica WRITE cai
	GLOBAL vol, TotalBuffer
	GLOBAL k1buf, k2buf, DCa
	RANGE cai0
	
	RANGE debug_flag
}

DEFINE NANN  4

UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(um) =	(micron)
	(mA) =	(milliamp)
	FARADAY =	(faraday)	(10000 coulomb)
	PI = (pi)	(1)
}

PARAMETER {
	DCa = 0.6			(um2/ms)
	
	
	k1buf	= 100			(/mM-ms)
	k2buf	= 0.1			(/ms)
	TotalBuffer = 0.003	(mM)
	cai0 = 50e-6 (mM)	
	

}

ASSIGNED {
	diam		(um)
	ica		(mA/cm2)
	cai		(mM)
	vol[NANN]	(1)	
	Kd		(/mM)
	B0		(mM)
}

STATE {
	ca[NANN]		(mM) <1e-6>	
	CaBuffer[NANN]	(mM)
	Buffer[NANN]	(mM)
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}

	cai = cai0
	Kd = k1buf/k2buf
	B0 = TotalBuffer/(1 + Kd*cai)

	FROM i=0 TO NANN-1 {
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
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
	COMPARTMENT i, diam*diam*vol[i] {ca CaBuffer Buffer}
	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vol[i] {ca}
	~ ca[0] << (-ica*PI*diam/(2*FARADAY))
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
	}
	cai = ca[0]
}