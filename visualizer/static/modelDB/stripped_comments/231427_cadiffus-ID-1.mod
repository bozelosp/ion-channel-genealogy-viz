NEURON {
	SUFFIX cadiffus
	USEION ca READ cai, ica WRITE cai
	GLOBAL vrat
}

DEFINE Nannuli  4

UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(um) =	(micron)
	(mA) =	(milliamp)
	FARADAY =	(faraday)	(10000 coulomb)
	PI = (pi)	(1)
}

PARAMETER {
	DCa = 	0.23		(um2/ms) 
}

ASSIGNED {
	diam	(um)
	ica		(mA/cm2)
	cai		(mM)
	vrat[Nannuli]		
                        
                        
	B0		(mM)
}

CONSTANT { volo = 1e10 (um2)}

STATE {
	ca[Nannuli]		(mM) <1e-6>	
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

	FROM i=0 TO Nannuli-1 {
		ca[i] = cai
	}
}

LOCAL frat[Nannuli]

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2			
	dr2 = r/(Nannuli-1)/2	
	vrat[0] = 0
	frat[0] = 2*r
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	
					
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2	
	}
}

LOCAL dsq, dsqvol	
			

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer}
	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}
	~ ca[0] << (-ica*PI*diam/(2*FARADAY))
	FROM i=0 TO Nannuli-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
cai = ca[0]
}