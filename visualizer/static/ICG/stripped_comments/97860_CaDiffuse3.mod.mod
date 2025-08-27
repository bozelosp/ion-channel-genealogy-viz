NEURON {
	SUFFIX CaDiffuse3
	USEION ca READ cai, ica WRITE cai
	GLOBAL vol, TotalBuffer, Pmax, beta, Dapp, mult
	RANGE cai0, d2cai, cai_prime, cal, d1cai_0, d1cai_1, counter
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
	Pmax = 2	(um/ms)
	beta = 0.001
	z = 2
	F = 9.64846e4	(coulomb)
	Dapp = .02	(um2/ms)
	mult = .000001
}

ASSIGNED {
	diam	(um)
	ica		(mA/cm2)
	icalcium	(mA/m2)
	dt		(ms)
	cai_prime	(mM/ms)
	cai		(mM)
	vol[NANN]	(1)	
	Kd		(/mM)
	B0		(mM)
	d1cai_0	(mM/um)
	d1cai_1	(mM/um)
	d2cai	(mM/um2)
	cal		(mM)
	counter
}

STATE {
	ca[NANN]		(mM) <1e-6>	
	CaBuffer[NANN]	(mM)
	Buffer[NANN]	(mM)
}

BREAKPOINT {
	ca[0]=cai
	SOLVE state METHOD sparse
}

LOCAL factors_done

INITIAL {
	counter = 0
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}
	
	Kd = k1buf/k2buf
	B0 = TotalBuffer/(1 + Kd*cai)

	FROM i=0 TO NANN-1 {
		ca[i] = cai*(1-i/NANN)
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
	}
	d2()
	rate()
	caistep()
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
	
	d2()
	rate()
	caistep()
}

PROCEDURE d2() {LOCAL radius, dr
	radius = diam/2
	dr = radius/(NANN-1)
	d1cai_0 = (ca[0] - ca[1])/dr
	d1cai_1 = (ca[1] - ca[2])/dr

	d2cai = (d1cai_0 - d1cai_1)/dr
}

PROCEDURE rate() {
	icalcium = ica * 1e4 (cm2/m2)
	cai_prime = (Dapp*d2cai) -((icalcium*4*beta)/(z*F*diam)) - ((Pmax*4*beta*ca[0])/diam)
}

PROCEDURE caistep() {
	counter = counter+1
	cai = ca[0] + cai_prime*dt
}