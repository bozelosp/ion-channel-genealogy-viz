NEURON {
	SUFFIX cadp
	USEION ca READ cao, cai, ica WRITE cai, ica
	RANGE ica_pump
	GLOBAL vol, TotalBuffer, TotalPump
	RANGE cai0
	THREADSAFE
}

DEFINE NANN  4

UNITS {
	(molar) =	(1/liter)
	(mol) = (1)
	(mM) =	(millimolar)
	(um) =	(micron)
	(mA) =	(milliamp)
	FARADAY =	(faraday)	(10000 coulomb)
	PI = (pi)	(1)
}

PARAMETER {
	DCa = 0.6		(um2/ms)
	
	
	k1buf	= 100			(/mM-ms)
	k2buf	= 0.1			(/ms)
	TotalBuffer = 0.003	(mM)
	cai0 = 50e-6 (mM)	
	
	k1 = 1 (/mM-ms)
	k2 = 0.005(/ms)
	k3 = 1 (/ms)
	k4 = 0.005 (/mM-ms)
	
	TotalPump = 1e-11 (mol/cm2)
	
}

ASSIGNED {
	diam		(um)
	ica		(mA/cm2)
	cai		(mM)
	vol[NANN]	(1)	
	Kd		(/mM)
	B0		(mM)
	
	cao (mM)
	ica_pmp (mA/cm2)
	parea (um)
}

CONSTANT{ volo = 1e10 (um2) }

STATE {
	ca[NANN]		(mM) <1e-6>	
	CaBuffer[NANN]	(mM)
	Buffer[NANN]	(mM)
	
	pump (mol/cm2)
	pumpca (mol/cm2)
}

BREAKPOINT {
	SOLVE state METHOD sparse
	ica = ica_pmp
}

LOCAL factors_done

INITIAL {
	MUTEXLOCK
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}
	MUTEXUNLOCK

	cai = cai0
	Kd = k1buf/k2buf
	B0 = TotalBuffer/(1 + Kd*cai)

	FROM i=0 TO NANN-1 {
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
	}
	
	parea = PI*diam
	pump = TotalPump/(1+(cai*k1/k2))
	pumpca = TotalPump - pump
	
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
	
	COMPARTMENT (1e10)*parea {pump pumpca}
	COMPARTMENT volo {cao}	
	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vol[i] {ca}
	
	
	~ ca[0] + pump <-> pumpca (k1*parea*(1e10), k2*parea*(1e10))
	~ pumpca <-> pump + cao  (k3*parea*(1e10), k4*parea*(1e10))
	CONSERVE pump + pumpca = TotalPump *parea * (1e10)
	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea
	
	
	
	~ ca[0] << (-(ica - ica_pmp)*PI*diam/(2*FARADAY))
	
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