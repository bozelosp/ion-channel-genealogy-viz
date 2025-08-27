NEURON {
  POINT_PROCESS Izhi2007b_dyn_thr
  RANGE C, k, vr, v0, u0, vt, vt0, vpeak, dyn_th, u, a, b, c, d, alpha, Vi, Vmin, Ka, Ki, tau, Iin, celltype, alive, cellid, verbose, derivtype, delta, t0
  NONSPECIFIC_CURRENT i
}


UNITS {
  (mV) = (millivolt)
  (uM) = (micrometer)
}


PARAMETER {
  C = 1 
  k = 0.7
  vr = -60 (mV)   
  v0 = -70 (mV)   
  u0 = 0          
  vt0 = -40 (mV)  
  vpeak = 35 (mV) 
  dyn_th = 0      
  a = 0.03
  b = -2
  c = -50
  d = 100
  Iin = 0
  alpha = 0.3
  Vi = -41.0
  Vmin = -37.0
  Ka = 7
  Ki = 8.75
  tau = 6.5
  celltype = 1 
  alive = 1 
  cellid = -1 
}


ASSIGNED {
  v (mV)
  i (nA)
  u (mV) 
  vt (mV)
  delta
  t0
  derivtype
}


INITIAL {
  u = u0
  if (dyn_th == 1) {
    vt = vt_thetainf(vr)
  } else {
    vt = vt0
  }
  derivtype=2
  net_send(0,1) 
}


BREAKPOINT {
  delta = t-t0 
  if (celltype<5) {
    if (t>0) {
      u = u + delta*a*(b*(v-vr)-u) 
      if (dyn_th == 1) {
        vt = vt + delta*(vt_thetainf(v) - vt)/tau
      }
    }
  }
  else {
     
     if (celltype==5) {
       if (v<d) { 
        u = u + delta*a*(0-u)
       }
       else { 
        u = u + delta*a*((0.025*(v-d)*(v-d)*(v-d))-u)
       }
     }

     
     if (celltype==6) {
       if (v>-65) {b=0}
       else {b=15}
       u = u + delta*a*(b*(v-vr)-u) 
     }
     
     
     if (celltype==7) {
       if (v>-65) {b=2}
       else {b=10}
       u = u + delta*a*(b*(v-vr)-u) 
     }
  }

  t0=t 
  i = -(k*(v-vr)*(v-vt) - u + Iin)/C/1000
}

FUNCTION vt_thetainf (v(mV)) (mV){
  vt_thetainf = alpha * (v - Vi) + Vmin + Ka * log ( 1 + exp((v - Vi)/Ki) )
}

FUNCTION derivfunc () {
  if (celltype==5 && derivtype==2) { 
    derivfunc = a*(0-u)
  } else if (celltype==5 && derivtype==1) { 
    derivfunc = a*((0.025*(v-d)*(v-d)*(v-d))-u)
  } else if (celltype==5) { 
    VERBATIM
    hoc_execerror("izhi2007b.mod ERRA: derivtype not set",0);
    ENDVERBATIM
  } else {
    derivfunc = a*(b*(v-vr)-u) 
  }
}


NET_RECEIVE (w) {
  
  if (flag == 1) { 
    if (celltype == 4) { 
      WATCH (v>(vpeak-0.1*u)) 2 
    } else if (celltype == 6) { 
      WATCH (v>(vpeak+0.1*u)) 2 
    } else { 
      WATCH (v>vpeak) 2 
    }
    
    if (celltype==6 || celltype==7) {
      WATCH (v> -65) 3 
      WATCH (v< -65) 4 
    }
    if (celltype==5) {
      WATCH (v> d) 3  
      WATCH (v< d) 4  
    }
    v = v0  
  
  } else if (flag == 2) { 
    if (alive) {net_event(t)} 
    
    if (celltype == 4) {
      v = c+0.04*u 
      if ((u+d)<670) {u=u+d} 
      else {u=670} 
     }  
    
    else if (celltype == 5) {
      v = c 
     }  
    
    else if (celltype == 6) {
      v = c-0.1*u 
      u = u+d 
     }  else {
      v = c 
      u = u+d 
     }
  
  } else if (flag == 3) { 
    
    if (celltype == 5)        { derivtype = 1 
    } else if (celltype == 6) { b=0
    } else if (celltype == 7) { b=2 
    }
  
  } else if (flag == 4) { 
    if (celltype == 5)        { derivtype = 2  
    } else if (celltype == 6) { b=15
    } else if (celltype == 7) { b=10
    }
  }
}