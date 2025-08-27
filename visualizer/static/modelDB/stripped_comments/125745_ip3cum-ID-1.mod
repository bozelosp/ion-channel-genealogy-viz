NEURON {
  SUFFIX ip3cum
  USEION ip3 READ ip3i, iip3 WRITE ip3i VALENCE 1
  GLOBAL vrat  
}

DEFINE Nannuli 4

UNITS {
  (molar) = (1/liter)
  (mM)    = (millimolar)
  (uM)    = (micromolar)
  (um)    = (micron)
  (mA)    = (milliamp)
  FARADAY = (faraday)  (coulomb)
  PI      = (pi)       (1)
}

PARAMETER {
  DIP3 = 0.283 (um2/ms)
  kdegr = 0.14e-3 (/ms)  
}

CONSTANT {
  ip3i0 = 0.16e-3 (mM)  
}

ASSIGNED {
  diam      (um)
  iip3      (mA/cm2)
  ip3i      (mM)
  vrat[Nannuli]  
                 
                 
}

STATE {
  
  
  ip3[Nannuli]       (mM) <1e-6>
}

BREAKPOINT { SOLVE state METHOD sparse }

LOCAL factors_done

INITIAL {
   if (factors_done == 0) {  
      factors_done = 1       
      factors()              
   }

  ip3i = ip3i0
  FROM i=0 TO Nannuli-1 {
    ip3[i] = ip3i
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
  COMPARTMENT i, diam*diam*vrat[i] {ip3 ip3i0}
  LONGITUDINAL_DIFFUSION i, DIP3*diam*diam*vrat[i] {ip3}
  ~ ip3[0] << (-iip3*PI*diam*(1e4)/FARADAY)  
  FROM i=0 TO Nannuli-2 {
    ~ ip3[i] <-> ip3[i+1]  (DIP3*frat[i+1], DIP3*frat[i+1])
  }
  dsq = diam*diam
  FROM i=0 TO Nannuli-1 {
    dsqvol = dsq*vrat[i]
    ~ ip3[i] <-> ip3i0  (kdegr*dsqvol, kdegr*dsqvol)
  }
  ip3i = ip3[0]
}