NEURON {
	SUFFIX cabalstore
	USEION ca READ cai, ica WRITE cai

	THREADSAFE 
	RANGE  cainit, fCa , icapump,icapumpmax,km, TotalBuffer, k1buf, k2buf, SCALE, shellfrac, DCa, vrat, frat, basal, cac, cast, tog
	POINTER ju_p, jcicr_p
}

DEFINE Nannuli 2 

UNITS {
	(molar) = (1/liter)
	(mM) =  (millimolar)
	(um) =  (micron)
	(mA) =  (milliamp)
	FARADAY = (faraday) (coulomb) 
	PI = (pi) (1)
}

PARAMETER {
        fCa = 0.05  (1)
        cainit = 0.0001 (mM)
        dt    (ms)
        celsius = 35  (degC)
        icapumpmax  = 0.00191  (mA/cm2) 
        km = 0.000500         (mM)
        DCa = 0.6 (um2/ms) 
        k1buf = 100 (/mM-ms) 
        k2buf = 0.1 (/ms)
        TotalBuffer = 0.03  (mM)
        basal = 0 (mM) 
        shellfrac = 0.25 
        SCALE = 1
        tog = 1

         }

ASSIGNED {
  diam  (um)
  ica   (mA/cm2)
  icapump (mA/cm2)
  vrat[Nannuli]
  Kd  (/mM)
  B0  (mM)
  cai (mM)
  cast (mM)
  cac (mM)
  jcicr_p
  ju_p
}

STATE {
  ca[Nannuli] (mM) <1e-10>
  CaBuffer[Nannuli] (mM) <1e-10>
  Buffer[Nannuli] (mM) <1e-10>
  castore[Nannuli] (mM) <1e-10>
}

BREAKPOINT {
	SOLVE states METHOD sparse
}

INITIAL{
  factors()
   
  cai=cainit
  
  Kd = k1buf
  B0 = TotalBuffer/(1+Kd*cainit)
  
  FROM i=0 TO Nannuli-1 {
    ca[i] = cai
    Buffer[i] = B0
    CaBuffer[i] = TotalBuffer - B0
    castore[i] = 0.2 
  }
  
  
}

LOCAL frat[Nannuli]
PROCEDURE factors() {
  LOCAL r, dr2, shsq
  if(Nannuli > 1) {
  r = 1/2 
  dr2 = r/(Nannuli-1)/2 
  
  vrat[0] = 0
  frat[0] = 2*r*PI
  FROM i=0 TO Nannuli-2 {
    vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2 
    r = r - dr2
    frat[i+1] = 2*PI*r/(2*dr2) 
                               
    r = r - dr2
    vrat[i+1] = PI*(r+dr2/2)*2*dr2 
    }
  } else { 
  frat[0] = 1
  vrat[0] = 1
  }
  if(Nannuli == 2){
     frat[1] = PI*(1.0-shellfrac)
     vrat[1] = 0.25*PI*(1.0-shellfrac)*(1.0-shellfrac)
     vrat[0] = 0.25*PI-vrat[1]
  
  }

  }
 
LOCAL dsq, dsqvol, casq



KINETIC states {
  COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer Buffer castore}
  icapump = icapumpmax*(1/(1 + km/ca[0])) 

  FROM i = 0 TO Nannuli-2 {
    ~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1]) 
  }
  dsq = diam*diam
  FROM i=0 TO Nannuli-1 {
    dsqvol = dsq*vrat[i]
    ~ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
  }
  dsqvol = dsq*vrat[1]
  ~ ca[Nannuli-1] <-> castore[Nannuli-1] (tog*ju_p*dsqvol,tog*jcicr_p*dsqvol)
  ~ ca[0] << (-SCALE*(ica +icapump)*PI*diam*(1e4)/(2*FARADAY))
  
  cai = ca[0]
  cac = ca[Nannuli-1]
  cast = castore[Nannuli-1]
}