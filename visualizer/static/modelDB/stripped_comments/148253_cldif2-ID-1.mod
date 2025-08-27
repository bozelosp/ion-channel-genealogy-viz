NEURON {
	SUFFIX cldifus2
	USEION cl READ icl WRITE cli VALENCE -1
	USEION hco3 READ hco3i, hco3o VALENCE -1
	GLOBAL vrat		
	RANGE cli0, vmax, leak, Kd
}

DEFINE Nannuli 4

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
}

PARAMETER {
	DCl = 2 (um2/ms) 
	leak = .000588235 (mM/ms) 	
	vmax = .005 (mM/ms)	
	Kd = 15 (mM)
	cli0 = 2 (mM)
	
	
	
	
}

ASSIGNED {
	diam 	(um)
	icl 	(mA/cm2)
	cli 	(mM)
	hco3i	(mM)
	hco3o	(mM)
	vrat[Nannuli]	
			
			
}

STATE {
	
	
	cl[Nannuli]	(mM) <1e-10>
}


BREAKPOINT { SOLVE state METHOD sparse}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {  	
		factors_done = 1	
		factors()		
	}
	cli = cli0
	FROM i=0 TO Nannuli-1 {
		cl[i] = cli
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

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {cl}
	LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}
	~ cl[0] << ((icl*PI*diam/FARADAY) + (leak - vmax*(cl[0]/(Kd + cl[0])))*diam*diam*vrat[0]) 
	FROM i=0 TO Nannuli-2 {
		~ cl[i] <-> cl[i+1]	(DCl*frat[i+1], DCl*frat[i+1])
	}
	cli = cl[0]
}