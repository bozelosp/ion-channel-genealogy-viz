INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
  (pA) = (picoamp)
  (mV) = (millivolt)
  (nS) = (nanosiemens)
}

NEURON {
  POINT_PROCESS AMPA16v8_noNC
  POINTER Glu
  RANGE LTP_ampaNbModFactor
  RANGE kass_re1   
  RANGE kdiss_re1  
  RANGE kass_re5
  RANGE kdiss_re5
  RANGE kass_re11
  RANGE kdiss_re11
  RANGE kass_re12
  RANGE kdiss_re12
  RANGE kass_re16
  RANGE kdiss_re16
  RANGE kass_re19
  RANGE kdiss_re19
  RANGE conduc_O2
  RANGE conduc_O3
  RANGE conduc_O4
  RANGE Erev_AMPA
  RANGE sumOpen
  RANGE PNa
  RANGE PK
  RANGE PCa
  RANGE ICa_AMPA
  RANGE INa_AMPA
  RANGE IK_AMPA
  RANGE nbAMPAR
  RANGE NewNbAMPAR
  RANGE Deact_factor
  RANGE Desens_factor
  RANGE kdiss_re16_Init
  RANGE kass_re11_Init
  RANGE kass_re12_Init
  RANGE position_AMPAR
  NONSPECIFIC_CURRENT i
  RANGE g
  RANGE v1
}

PARAMETER {
  kass_re1 = 10.0     
  kdiss_re1 = 7.0     
  kass_re5 = 10.0     
  kdiss_re5 = 0.00041 
  kdiss_re11 = 0.001  
  kdiss_re12 = 0.017  
  kass_re16 = 0.55    
  kass_re19 = 0.2     
  kdiss_re19 = 0.035  
  conduc_O2 = 9.0
  conduc_O3 = 15.0
  conduc_O4 = 21.0
  Erev_AMPA = 0.0
  PNa = 50.0
  PK = 49.5
  PCa = 0.5
  nbAMPAR = 80
  Deact_factor = 1.0
  Desens_factor = 1.0
  kdiss_re16_Init = 0.3    
  kass_re11_Init = 3.3e-06 
  kass_re12_Init = 0.42    
  position_AMPAR = 60.0
  LTP_ampaNbModFactor = 1
  v
  Glu
  v1
}

STATE {
  R0
  R1
  R2
  R3
  R4
  D0
  D1
  D2
  D3
  D4
  E2
  E3
  E4
  O2
  O3
  O4
}

INITIAL {
  R0 = 1.0
  R1 = 0.0
  R2 = 0.0
  R3 = 0.0
  R4 = 0.0
  D0 = 0.0
  D1 = 0.0
  D2 = 0.0
  D3 = 0.0
  D4 = 0.0
  E2 = 0.0
  E3 = 0.0
  E4 = 0.0
  O2 = 0.0
  O3 = 0.0
  O4 = 0.0
}
ASSIGNED{
  kdiss_re16
  kass_re11
  kass_re12
  NewNbAMPAR
  sumOpen
  i
  INa_AMPA
  IK_AMPA
  ICa_AMPA
  g

}


BREAKPOINT {
  SOLVE states METHOD derivimplicit

  kdiss_re16 = kdiss_re16_Init / Deact_factor
  kass_re11 = kass_re11_Init / Desens_factor
  kass_re12 = kass_re12_Init / Desens_factor

  


  NewNbAMPAR = nbAMPAR * (16/40) 

  sumOpen = O2 + O3 + O4
  g= (conduc_O2 * O2 + conduc_O3 * O3 + conduc_O4 * O4) * NewNbAMPAR * 1e-3 
  i= (conduc_O2 * O2 + conduc_O3 * O3 + conduc_O4 * O4) * (v- Erev_AMPA) * NewNbAMPAR * 1e-3 

  INa_AMPA = PNa / 100 * i
  IK_AMPA = PK / 100 * i
  ICa_AMPA = PCa / 100 * i

  v1 = v
}

DERIVATIVE states {
  LOCAL dummy ,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p1a,p1b,p1c,p1d,p1e,p1f,p20
  p1 = kass_re1
  p2 = kdiss_re1
  p3 = kass_re5
  p4 = kdiss_re5
  p5 = kass_re11
  p6 = kdiss_re11
  p7 = kass_re12
  p8 = kdiss_re12
  p9 = kass_re16
  pa = kdiss_re16
  pb = kass_re19
  pc = kdiss_re19
  pd = conduc_O2
  pe = conduc_O3
  pf = conduc_O4
  p10 = Erev_AMPA


  p13 = PNa
  p14 = PK
  p15 = PCa
  p16 = ICa_AMPA
  p17 = INa_AMPA
  p18 = IK_AMPA
  p19 = nbAMPAR
  p1a = NewNbAMPAR
  p1b = Deact_factor
  p1c = Desens_factor
  p1d = kdiss_re16_Init
  p1e = kass_re11_Init
  p1f = kass_re12_Init
  p20 = position_AMPAR

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  R0' =  - (4*p1*R0*Glu - 1*p2*R1) - (4*p5*R0 - p6*D0) 
  R1' =  (4*p1*R0*Glu - 1*p2*R1) - (3*p1*R1*Glu - 2*p2*R2) - (1*p1f*R1 - p8*D1) 
  R2' =  (3*p1*R1*Glu - 2*p2*R2) - (2*p1*R2*Glu - 3*p2*R3) - (2*p1f*R2 - p8*D2) - (2*p9*R2 - pa*O2) 
  R3' =  (2*p1*R2*Glu - 3*p2*R3) - (1*p1*R3*Glu - 4*p2*R4) - (3*p1f*R3 - p8*D3) - (3*p9*R3 - pa*O3) 
  R4' =  (1*p1*R3*Glu - 4*p2*R4) - (4*p1f*R4 - p8*D4) - (4*p9*R4 - pa*O4) 
  D0' =  - (3*p3*D0*Glu - p4*D1) + (4*p5*R0 - p6*D0) 
  D1' =  (3*p3*D0*Glu - p4*D1) - (3*p1*D1*Glu - p2*D2) + (1*p1f*R1 - p8*D1) 
  D2' =  (3*p1*D1*Glu - p2*D2) - (2*p1*D2*Glu - 2*p2*D3) + (2*p1f*R2 - p8*D2) - (1*pb*D2 - pc*E2) 
  D3' =  (2*p1*D2*Glu - 2*p2*D3) - (1*p1*D3*Glu - 3*p2*D4) + (3*p1f*R3 - p8*D3) - (2*pb*D3 - pc*E3) 
  D4' =  (1*p1*D3*Glu - 3*p2*D4) + (4*p1f*R4 - p8*D4) - (3*pb*D4 - pc*E4) 
  E2' =  - (2*p1*E2*Glu - p2*E3) + (1*pb*D2 - pc*E2) 
  E3' =  (2*p1*E2*Glu - p2*E3) - (p1*E3*Glu - 2*p2*E4) + (2*pb*D3 - pc*E3) 
  E4' =  (p1*E3*Glu - 2*p2*E4) + (3*pb*D4 - pc*E4) 
  O2' =  (2*p9*R2 - pa*O2) 
  O3' =  (3*p9*R3 - pa*O3) 
  O4' =  (4*p9*R4 - pa*O4) 
}