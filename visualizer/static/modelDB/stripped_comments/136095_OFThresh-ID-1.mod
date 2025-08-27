UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
    (uS) = (microsiemens)
}

NEURON {
  POINT_PROCESS OFTH
  USEION other WRITE iother VALENCE 1.0

  RANGE gkbase             
  RANGE gkmin              
  RANGE vth                
  RANGE tauvtha, vthinc    
  RANGE taugka, gkinc      
  RANGE ik,ek              
  RANGE tauk               
  RANGE i,spkht            
  RANGE refrac             
  RANGE inrefrac           
  RANGE apdur              
  RANGE gna,ena,ina,gnamax 

  RANGE verbose
  GLOBAL checkref          
}

ASSIGNED { 
  v (mV)
  iother (nA)
  inrefrac
}

STATE { gk vthadapt gkadapt }

PARAMETER {
  gkbase=0.060(uS) 
  taugka=100 (ms) 
  kadapt=0.007(uS) 
  spkht = 55(mV)
  tauk=2.3 (ms) 
  ek = -70(mV)
  vth = -40(mV)
  refrac = 2.7(ms)
  verbose = 0
  apdur = 0.9 (ms)
  gkinc = 0.006(uS)
  tauvtha = 1(ms)
  vthinc = 0
  gna = 0(uS)
  ena = 55(mV)
  ina = 0(nA)
  ik = 0(nA)
  i = 0(nA)
  gnamax = .300(uS)
  gkmin = 0.00001(uS)
  checkref = 1
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  if( gk < gkmin ) { gk = gkmin }
  if( gkadapt < gkbase ) { gkadapt = gkbase }
  if( vthadapt < vth ) { vthadapt = vth }
  iassign()
}

INITIAL {
  net_send(0,1)
  gk = 0(uS)
  gkadapt = gkbase
  vthadapt = vth
  gna = 0(uS)
  ina = 0
  ik = 0  
  i = 0
  iother = 0
  inrefrac = 0
}

DERIVATIVE states {  
  gk' = -gk/tauk
  gkadapt' = (gkbase - gkadapt)/taugka
  vthadapt' = (vth - vthadapt)/tauvtha
}

PROCEDURE iassign () {
  ik = gk*(v-ek)
  ina = gna*(v-ena)
  i = ik + ina
  iother = i
}

NET_RECEIVE (w) {
  if (flag==1) {
    WATCH (v > vthadapt) 2
  } else if (flag==2 && !inrefrac) {  
      net_event(t)           
      net_send(apdur,3)      
      net_send(refrac,4)     
      inrefrac=1             
      gkadapt = gkadapt + gkinc 
      vthadapt = vthadapt + vthinc 
      gna = gnamax 
      if (verbose) { printf("spike at t=%g\n",t) }
  } else if(flag==3) {   
    gk = gkadapt           
    gna = 0 
    if (verbose) { printf("end of action potential @ t = %g\n",t) }
  } else if(flag==4) {   
    inrefrac = 0           
    if (verbose) { printf("refrac over @ t = %g\n",t) }
    
    if(checkref && v > vthadapt) { net_send(0,2) }
  } else if (flag==0 && w>0) {
    net_event(t)           
  } else if (flag==2 && inrefrac && verbose ) { printf("in refrac @ t = %g, no spike\n",t) }
}

FUNCTION fflag () { fflag=1 }

PROCEDURE version () {
  printf("$Id
}