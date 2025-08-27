NEURON {
	SUFFIX CaDiffuse2
	USEION ca READ cai, ica 
	POINTER cal
	GLOBAL vol, TotalBuffer
	RANGE cai0, d2cai
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
	diam	(um)
	ica		(mA/cm2)
	cai		(mM)
	vol[NANN]	(1)	
	Kd		(/mM)
	B0		(mM)
	d1cai_0	(mM/um)
	d1cai_1	(mM/um)
	d2cai	(mM/um2)
	cal		(mM)
}

STATE {
	ca[NANN]		(mM) <1e-6>	
	CaBuffer[NANN]	(mM)
	Buffer[NANN]	(mM)
}

BREAKPOINT {
	ca[0]=cal
	SOLVE state METHOD sparse
	d2()
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}

	Kd = k1buf/k2buf
	B0 = TotalBuffer/(1 + Kd*cal)

	FROM i=0 TO NANN-1 {
		ca[i] = cal*(1-i/NANN)
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
	}
	d2()
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
	
	
	
	
	cal = ca[0]
}

PROCEDURE d2() {LOCAL radius, dr
	radius = diam/2
	dr = radius/(NANN-1)
	d1cai_0 = (ca[0] - ca[1])/dr
	d1cai_1 = (ca[1] - ca[2])/dr

	d2cai = (d1cai_0 - d1cai_1)/dr
}