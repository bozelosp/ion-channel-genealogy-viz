NEURON {
  SUFFIX caquant
  USEION ca READ ica, cai
  
  
  GLOBAL mode
  RANGE iinit, imax, timax
  RANGE ihalf, t0i, t1i, hwi
  RANGE cinit, cmax, tcmax

  RANGE svr 
  RANGE qapprox, cmaxp

  RANGE vinit, vmax, tvmax
  RANGE vhalf, t0v, t1v, hwv
  RANGE vxt
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mM) = (milli/liter)
  F = (faraday) (coulombs)
}

PARAMETER {
  
  
  mode = 1 (1) 
}

ASSIGNED {
  ica (mA/cm2)
  iinit (mA/cm2) 
  imax (mA/cm2)  
  timax (ms)  
  ihalf (mA/cm2) 
  t0i (ms)    
  t1i (ms)    
  hwi (ms)    
  findwhichi (1) 

  cai (mM)
  cinit (mM) 
  cmax (mM)  
  tcmax (ms)  

  svr (/micron)  
            

  area (micron2) 
  qapprox (picocoulomb)
  cmaxp (mM)

  v (mV)
  vinit (mV) 
  vmax (mV)  
  tvmax (ms)  
  vhalf (mV) 
  t0v (ms)    
  t1v (ms)    
  hwv (ms)    
  findwhichv (1) 
  vxt (ms mV) 
}

INITIAL {
  if (mode==1) { 

    iinit = ica

    imax = ica
    timax = t
    ihalf = ica
    t0i = -1 (ms) 
    t1i = -1 (ms)
    hwi = -1 (ms)

    cinit = cai
    cmax = cai
    tcmax = t

    vinit = v
    vmax = v
    tvmax = t
    vhalf = v
    t0v = -1 (ms) 
    t1v = -1 (ms)
    hwv = -1 (ms)

  } else if (mode==2) { 

    ihalf = (iinit + imax)/2

    findwhichi = 0 
    vhalf = (vinit + vmax)/2
    findwhichv = 0 
  }
}


PROCEDURE findix() {
  if (findwhichi==0) {
    if (ica < ihalf) {
      t0i = t
      findwhichi = 1
    }
  } else if (findwhichi==1) {
    if (ica > ihalf) { 
      t1i = t
      hwi = t1i-t0i
      findwhichi = 2 
    }
  }
}


PROCEDURE findvx() {
  if (findwhichv == 0) {
    if (v > vhalf) {
      t0v = t
      findwhichv = 1
    }
  } else if (findwhichv == 1) {
    if (v < vhalf) { 
      t1v = t
      hwv = t1v-t0v
      findwhichv = 2 
    }
  }
}







BREAKPOINT { }


BEFORE STEP { 
  if (mode==1) { 
    if (ica<imax) {
      imax = ica
      timax = t
    }
    if (cai>cmax) {
      cmax = cai
      tcmax = t
    }
    if (v>vmax) {
      vmax = v
      tvmax = t
    }
  } else if (mode==2) {
    findix()
    qapprox = -(0.01)*area*(imax-iinit)*hwi
    cmaxp = cinit - (10000)*svr*(imax-iinit)*hwi/2/F
    findvx()
    vxt = (vmax - vinit)*hwv
  }
}