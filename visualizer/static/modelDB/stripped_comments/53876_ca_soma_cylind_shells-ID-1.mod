NEURON {
SUFFIX ca_soma_cyl
USEION ca READ cai, ica WRITE cai
GLOBAL vrat 

RANGE K2f_ex, K2f_ATPase, B
RANGE i_Na_Ca_ex, i_ATPase, I	
				
				
				
}
DEFINE Nannuli 50 
UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
}

PARAMETER {
	DCa = 0.6 (um2/ms)
	
	k1buf = 30  (/mM ms)  
	k2buf = 0.03 (/ms) 
		   
	B = 0.025 (mM) 
	
	
	
	

	K2f_ex = 0 
	
			
	cao = 2 (mM) 
		   
	E_1  = 0.01315 (/mV) 
	E_2 = 0.0255 (/mV)   
	nai = 7.6 (mM)	     
	nao = 152 (mM)	     
	K2f_ATPase = 0 
	f_ATPase = 100 (/mM ms)	
	b_ATPase = 0.005 (/ms)	
}

ASSIGNED {
	v (mV)
	diam (um)
	ica (mA/cm2)
	I (mA/cm2)	
	i_Na_Ca_ex (mA/cm2)
	i_ATPase (mA/cm2)
	cai (mM)
	vrat[Nannuli] (1) 
	
	
	
	Kd (/mM)
	B0 (mM)
}
STATE {
	
	
	ca[Nannuli] (mM) <1e10>
	CaBuffer[Nannuli] (mM)
	Buffer[Nannuli] (mM)
	n (1)
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	i_Na_Ca_ex = -K2f_ex * (nai^3 * cao * exp(E_1 * v) - nao^3 * cai * exp(-E_2*v))
	i_ATPase = K2f_ATPase * n
	I= ica	
	SOLVE state METHOD sparse
}

DERIVATIVE states {
	
	n' = f_ATPase * cai * (1 - n) - b_ATPase * n
}

LOCAL factors_done
INITIAL {
	if (factors_done == 0) { 
		factors_done = 1 
		factors() 
	}

	n = f_ATPase * cai / (f_ATPase * cai + b_ATPase)

	Kd = k1buf/k2buf
	B0 = B/(1 + Kd*cai)
	FROM i=0 TO Nannuli-1 {
		cai = 5e-5	
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = B - B0
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
		vrat[i] = vrat[i] + PI*(r-dr2/2)*4*dr2 
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2) 
		
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2 
	}
}
LOCAL dsq, dsqvol 

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer}
	
	~ ca[0] << ( (-ica - i_Na_Ca_ex - i_ATPase)*PI*diam/(2*FARADAY)) 
			
			
			
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