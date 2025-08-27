? interface
NEURON {
	SUFFIX cab
	NONSPECIFIC_CURRENT idummy 
	GLOBAL vol, frat
	RANGE in
}





DEFINE NANN  99

UNITS {
	(mV)	= (millivolt)
	(molar) = (1/liter)
	(uM)	= (micromolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(mol)	= (1)
	(pmol)	= (picomol)
	FARADAY = (faraday)	 (coulomb)
	PI	= (pi)		(1)
	R 	= (k-mole)	(joule/degC)
}

PARAMETER {
	DFree = .6	(um2/ms) <0,1e9>
	beta = 50		<0, 1e9>

	ke = 1e-3	(cm/s)	

	dur = 1		(ms) 
	t2 = 1e9	(ms)
	amp = 1000	(pmol/cm2-s) 

	dx = 100	(angstrom) <1, 1e9> 
	ca0 = 0		(um)
}

ASSIGNED {
	celsius		(degC)
	diam		(um)
	vol[NANN]	(um2)	
	frat[NANN] 	()	
	in		(pmol/cm2-s)
	idummy (mA/cm2)
}


STATE {
	ca[NANN]	(uM) <1e-6>
}

INITIAL {LOCAL total
	parms()
	FROM i=0 TO NANN-1 {
		ca[i] = ca0
	}
}

BREAKPOINT {
	SOLVE state METHOD sparse
	if (at_time(0)) {
		in = amp
	}
	if (at_time(dur)) {
		in = 0
	}
	if (at_time(t2)) {
		in = amp
	}
	if (at_time(t2+dur)) {
		in = 0
	}
	idummy = 0 
}


PROCEDURE coord() {
	LOCAL r, dr2
	
	
	
	
	
	
	r = diam/2			
	dr2 = dx/2*(1e-4)		
	vol[0] = 0
	frat[0] = 0 
	FROM i=0 TO NANN-2 {
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2	
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	
					
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2	
	}
	
	vol[NANN-1] = vol[NANN-1] + PI*r*r
}

KINETIC state {
	COMPARTMENT i, (1+beta)*vol[i]*1(um) {ca}
? kinetics
	~ ca[0] -> (ke*PI*diam*(10)*1(um)) 
	~ ca[0] << (in*PI*diam*(1e-2)*1(um)) 
	
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DFree*frat[i+1]*1(um), DFree*frat[i+1]*1(um))
	}
}
	
PROCEDURE parms() {
	coord()
}