NEURON {
SUFFIX ca_soma
USEION ca READ cai, ica WRITE cai
NONSPECIFIC_CURRENT ifake 
GLOBAL vrat 

GLOBAL factor	
RANGE K2f_ex, K2f_ATPase, B	
RANGE i_Na_Ca_ex, i_ATPase, I	
				
				
				
				
}
DEFINE Nthin	5	
			
			
			
			
DEFINE Nthick	21	

DEFINE Nshells	25	

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
	
	
	
	

	K2f_ex =  6.3680304e-08 (mA/cm2 mM4) 
	
			
	cao = 2 (mM) 
		   
	E_1  = 0.01315 (/mV) 
	E_2 = 0.0255 (/mV)   
	nai = 7.6 (mM)	     
	nao = 152 (mM)	     
	K2f_ATPase =0.00051137214  (mA/cm2) 
							 
	f_ATPase = 100 (/mM ms)	
	b_ATPase = 0.005 (/ms)	
	factor = 1 (1)	
	mM2M = 1e-3 (1)	
}

ASSIGNED {
	v (mV)
	diam (um)
	ica (mA/cm2)
	I (mA/cm2)	
	i_Na_Ca_ex (mA/cm2)
	i_ATPase (mA/cm2)
	cai (mM)
	vrat[Nshells] (um3)	
	Kd (/mM)
	B0 (mM)
	ifake (mA/cm2) 
		       
}
STATE {
	
	
	ca[Nshells] (mM) <1e-6>
	CaBuffer[Nshells] (mM) <1e-6>
	Buffer[Nshells] (mM) <1e-6>
	n (1)
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	SOLVE state METHOD sparse
	i_Na_Ca_ex = -K2f_ex * (nai^3 * cao * exp(E_1 * v) - nao^3 * cai * exp(-E_2*v))
	i_ATPase = K2f_ATPase * n
	I= ica	
	ifake=0 

}

DERIVATIVE states {
	
	n' = f_ATPase * cai * (1 - n) - b_ATPase * n
}

INITIAL {
	factors() 
	n = f_ATPase * cai / (f_ATPase * cai + b_ATPase)

	Kd = k1buf/k2buf
	B0 = B/(1 + Kd*cai)
	FROM i=0 TO Nshells-1 {
		cai = 5e-5 
		ca[i] = cai
		Buffer[i] = B0
		CaBuffer[i] = B - B0
	}
}
LOCAL frat[Nshells] 
LOCAL radii[Nthick] 
LOCAL drthick[Nthick] 
	
	
LOCAL drthick2[Nthick]	
LOCAL i_	

PROCEDURE factors() {
	LOCAL r, drthin, drthin2	
	r = diam/2 
	radii[0] = r
	FROM i=1 TO Nthick-1 {
		radii[i] = ( radii[i-1]^3 - (1/21)*r^3 )^(1/3)
		drthick[i-1] = radii[i-1] - radii[i]
		drthick2[i-1]=drthick[i-1]/2
	}




	drthin = 2 * (drthick[0]) / ( 2*Nthin - 1 )	
	drthin2 = drthin/2 
	
	vrat[0] = 0
	frat[0] = 2*r	
	
	FROM i=0 TO Nthin-2 {
		
		vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthin2-3*r*drthin2^2+drthin2^3)
		r = r - drthin2
		frat[i+1] = 4*PI*r*r/drthin
		r = r-drthin2
		vrat[i+1] = (4*PI/3) * (3*r^2*drthin2 + 3*r*drthin2^2 + drthin2^3)
	}
	
	FROM i=Nthin-1 TO Nshells-2 {
		i_ = 1	
			
		if (i==Nthin-1) {
			
			
			vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthin2-3*r*drthin2^2+drthin2^3)
			r = r - drthin2
			frat[i+1] = 4*PI*r*r/(drthin2+drthick2[i_])
			r = r - drthick2[i_]
			vrat[i+1] = (4*PI/3) * (3*r^2*drthick2[i_] + 3*r*drthick2[i_]^2 + drthick2[i_]^3)
		} else {
			vrat[i] = vrat[i] + (4*PI/3) * (3*r^2*drthick2[i_]-3*r*drthick2[i_]^2+drthick2[i_]^3)
			r = r - drthick2[i_]
			frat[i+1] = 4*PI*r*r/drthick[i_]
			r = r-drthick2[i_]
			vrat[i+1] = (4*PI/3) * (3*r^2*drthick2[i_] + 3*r*drthick2[i_]^2 + drthick2[i_]^3)
		}
		i_ = i_ + 1
	}
}



KINETIC state {
	COMPARTMENT i, vrat[i] {ca CaBuffer Buffer}
	
	~ ca[0] << ((-ica - i_Na_Ca_ex - i_ATPase)*PI*diam*diam /(2*FARADAY)) 
			
			
			
	FROM i=0 TO Nshells-2 {
		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
	}

	FROM i=0 TO Nshells-1 {

		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*vrat[i], k2buf*vrat[i])
	}
	cai = ca[0]
}