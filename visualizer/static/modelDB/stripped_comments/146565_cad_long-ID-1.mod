NEURON {
	SUFFIX cad
	USEION ca READ cao, cai, ica WRITE cai, ica
	RANGE ica_pmp
	GLOBAL vrat, TotalBuffer, TotalPump

}
DEFINE Nannuli 4

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
 	(nA) 	= (nanoamp)

}
CONSTANT { volo = 1e10 (um2) }
PARAMETER {
	DCa = 0.6 (um2/ms)
	k1buf = 100 (/mM-ms) 
	k2buf = 0.1 (/ms)
	TotalBuffer = 0.003 (mM)
	k1 = 1 (/mM-ms)
	k2 = 0.005 (/ms)
	k3 = 1 (/ms)
	k4 = 0.005 (/mM-ms)
	
	TotalPump = 1e-13(mol/cm2)
}
ASSIGNED {
	diam (um)
	
	ica (mA/cm2)
	cai (mM)
	vrat[Nannuli] 
	
	
	Kd (/mM)
	B0 (mM)
	cao (mM)
	ica_pmp (mA/cm2)
	parea (um)
}
STATE {
	
	
	ca[Nannuli] (mM) <1e-10>
	CaBuffer[Nannuli] (mM)
	Buffer[Nannuli] (mM)
	pump (mol/cm2)
	pumpca (mol/cm2)
}
BREAKPOINT {
	SOLVE state METHOD sparse
	ica = ica_pmp
}
LOCAL factors_done
INITIAL {
	parea = PI*diam
	if (factors_done == 0) {
		factors_done = 1
		factors()
	}
	Kd = k1buf/k2buf
	B0 = TotalBuffer/(1 + Kd*cai)
	FROM i=0 TO Nannuli-1 {
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = TotalBuffer - B0
	}
	pump = TotalPump/(1 + (cai*k1/k2))
	pumpca = TotalPump - pump
	ica=0
	cai=1e-5
	cao=2
	ica_pmp=0
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
	COMPARTMENT (1e10)*parea {pump pumpca}
	COMPARTMENT volo {cao}
	LONGITUDINAL_DIFFUSION i, DCa*vrat[i] {ca}
	
	~ ca[0] + pump <-> pumpca (k1*parea*(1e10), k2*parea*(1e10))
	~ pumpca <-> pump + cao (k3*parea*(1e10), k4*parea*(1e10))
	CONSERVE pump + pumpca = TotalPump * parea * (1e10)
	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

	
	~ ca[0] << (-(ica)*PI*diam/(2*FARADAY)) 
	FROM i=0 TO Nannuli-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO Nannuli-1 {
		dsqvol = dsq*vrat[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
	}
	cai = ca[0]
}