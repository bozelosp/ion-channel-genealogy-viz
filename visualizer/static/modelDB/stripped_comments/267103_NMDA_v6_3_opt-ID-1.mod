UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
}

NEURON {
  POINT_PROCESS NMDA_v6_3_opt
  POINTER Glu
  RANGE open_total
  RANGE cond
  RANGE g1
  RANGE g2
  RANGE alpha
  RANGE temp
  RANGE Vrev
  RANGE open_Mg
  NONSPECIFIC_CURRENT i
  RANGE nbNMDAR
  RANGE current_Ca
  RANGE perm_Ca
  RANGE current_K
  RANGE perm_K
  RANGE current_Na
  RANGE perm_Na
  RANGE ke
  RANGE k_e
  RANGE kg
  RANGE k_g
  RANGE gb
  RANGE g_b
  RANGE ga
  RANGE g_a
  RANGE don
  RANGE doff
  RANGE beta2
  RANGE alpha2
  RANGE beta1
  RANGE alpha1
  RANGE g
  RANGE v1
}

PARAMETER {

  cond = 0.0
  g1 = 40.0
  g2 = 247.0
  alpha = 0.01
  temp = 0.0
  
  Vrev = 0  
  open_Mg = 0.0
  nbNMDAR = 1.0
  current_Ca = 0.0
  perm_Ca = 0.1
  current_K = 0.0
  perm_K = 0.45
  current_Na = 0.0
  perm_Na = 0.45
  ke = 1.0
  k_e = 0.0263
  kg = 10.0
  k_g = 0.0291
  gb = 0.168153226487
  g_b = 0.263214630793
  ga = 0.1
  g_a = 217.59688772
  don = 0.0421976169622
  doff = 0.0128627223774
  beta2 = 7.11918962116
  alpha2 = 4.36739898423
  beta1 = 3.5
  alpha1 = 0.174379055622
  Gly = 0.02
  Mg = 1
  Glu
  v
  v1
  g
}

STATE {
  R
  R_Glu
  R_2Glu
  R_Gly
  R_2Gly
  R_Glu_Gly
  R_Glu_2Gly
  R_2Glu_Gly
  R_2Glu_2Gly
  State6
  State5
  State4
  Desensitized
  Open1
  Open2
}

INITIAL {
  R = 1.0
  R_Glu = 0.0
  R_2Glu = 0.0
  R_Gly = 0.0
  R_2Gly = 0.0
  R_Glu_Gly = 0.0
  R_Glu_2Gly = 0.0
  R_2Glu_Gly = 0.0
  R_2Glu_2Gly = 0.0
  State6 = 0.0
  State5 = 0.0
  State4 = 0.0
  Desensitized = 0.0
  Open1 = 0.0
  Open2 = 0.0
}

ASSIGNED{

  i (nA)
  open_total
}

BREAKPOINT {
  SOLVE states METHOD derivimplicit
  

  
  open_total = Open2 + Open1

  
  cond = g1 + ((g2 - g1) / (1 + exp(alpha * v * 0.001(1/mV))))

  
  open_Mg = (Open1 + Open2) / (1 + exp(-62 * v * 0.001(1/mV)) * Mg / 3.57)
  
  
  temp = cond * (v - Vrev) * 0.001 * open_Mg

  
  i = temp * nbNMDAR

  
  current_Ca = perm_Ca * i

  
  current_K = perm_K * i

  
  current_Na = perm_Na * i
  
  g = cond * open_Mg * 0.001
  v1 = v
}

DERIVATIVE states {
  LOCAL dummy ,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p1a,p1b,p1c,p1d,p1e
  p1 = open_total
  p2 = cond
  p3 = g1
  p4 = g2
  p5 = alpha
  p6 = temp
  p7 = Vrev
  p8 = open_Mg
  pa = nbNMDAR
  pb = current_Ca
  pc = perm_Ca
  pd = current_K
  pe = perm_K
  pf = current_Na
  p10 = perm_Na
  p11 = ke
  p12 = k_e
  p13 = kg
  p14 = k_g
  p15 = gb
  p16 = g_b
  p17 = ga
  p18 = g_a
  p19 = don
  p1a = doff
  p1b = beta2
  p1c = alpha2
  p1d = beta1
  p1e = alpha1

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  R' =  - (2*p11*R*Glu - p12*R_Glu) - (2*p13*R*Gly - p14*R_Gly) 
  R_Glu' =  (2*p11*R*Glu - p12*R_Glu) - (p11*R_Glu*Glu - 2*p12*R_2Glu) - (2*p13*R_Glu*Gly - p14*R_Glu_Gly) 
  R_2Glu' =  (p11*R_Glu*Glu - 2*p12*R_2Glu) - (2*p13*R_2Glu*Gly - p14*R_2Glu_Gly) 
  
  R_Gly' =  (2*p13*R*Gly - p14*R_Gly) - (p13*R_Gly*Gly - 2*p14*R_2Gly) - (2*p11*R_Gly*Glu - p12*R_Glu_Gly) 
  R_2Gly' =  (p13*R_Gly*Gly - 2*p14*R_2Gly) - (2*p11*R_2Gly*Glu - p12*R_Glu_2Gly) 
  
  R_Glu_Gly' =  (2*p13*R_Glu*Gly - p14*R_Glu_Gly) - (p13*R_Glu_Gly*Gly - 2*p14*R_Glu_2Gly) + (2*p11*R_Gly*Glu - p12*R_Glu_Gly) - (p11*R_Glu_Gly*Glu - 2*p12*R_2Glu_Gly) 
  R_Glu_2Gly' =  (p13*R_Glu_Gly*Gly - 2*p14*R_Glu_2Gly) + (2*p11*R_2Gly*Glu - p12*R_Glu_2Gly) - (p11*R_Glu_2Gly*Glu - 2*p12*R_2Glu_2Gly) 
  R_2Glu_Gly' =  (2*p13*R_2Glu*Gly - p14*R_2Glu_Gly) - (p13*R_2Glu_Gly*Gly - 2*p14*R_2Glu_2Gly) + (p11*R_Glu_Gly*Glu - 2*p12*R_2Glu_Gly) 
  R_2Glu_2Gly' =  (p13*R_2Glu_Gly*Gly - 2*p14*R_2Glu_2Gly) + (p11*R_Glu_2Gly*Glu - 2*p12*R_2Glu_2Gly) - (p15*R_2Glu_2Gly - p16*State6) - (p17*R_2Glu_2Gly - p18*State5) 
  State6' =  (p15*R_2Glu_2Gly - p16*State6) - (p17*State6 - p18*State4) - (p19*State6 - p1a*Desensitized) 
  State5' =  (p17*R_2Glu_2Gly - p18*State5) - (p15*State5 - p16*State4) - (p1b*State5 - p1c*Open2) 
  State4' =  (p15*State5 - p16*State4) + (p17*State6 - p18*State4) - (p1d*State4 - p1e*Open1) 
  Desensitized' =  (p19*State6 - p1a*Desensitized) 
  Open1' =  (p1d*State4 - p1e*Open1) 
  Open2' =  (p1b*State5 - p1c*Open2) 
}