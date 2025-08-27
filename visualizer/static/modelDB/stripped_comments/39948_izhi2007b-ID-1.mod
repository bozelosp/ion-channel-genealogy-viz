NEURON {
  POINT_PROCESS Izhi2007b
  RANGE C, k, vr, vt, vpeak, u, a, b, c, d, Iin, celltype, alive, cellid, verbose, derivtype, delta, t0
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
  vt = -40 (mV) 
  vpeak = 35 (mV) 
  a = 0.03
  b = -2
  c = -50
  d = 100
  Iin = 0
  celltype = 1 
  alive = 1 
  cellid = -1 
}


ASSIGNED {
  v (mV)
  i (nA)
  u (mV) 
  delta
  t0
  derivtype
}


INITIAL {
  u = 0.0
  derivtype=2
  net_send(0,1) 
}


BREAKPOINT {
  delta = t-t0 
  if (celltype<5) {
    u = u + delta*a*(b*(v-vr)-u) 
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
    v = vr  
  
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