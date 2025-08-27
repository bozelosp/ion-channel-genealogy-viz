NEURON {
    POINT_PROCESS IzhiGS
    POINTER vpreGraded   
  RANGE C, k, vr, vt, vpeak, a, b, c, d, Iin, tauAMPA, tauGraded, celltype, alive, cellid, verbose
  RANGE V, u, gAMPA, gGraded, gbarGraded, I, eAMPA, eGraded, vmidGraded, vslopeGraded, sinfGraded
  RANGE factor, eventflag, delta, t0
}


UNITS {
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (um) = (micrometer)
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
  Vpre = 0
  gbarGraded = 0  (uS)
  tauAMPA = 5 (ms) 
  tauGraded = 4.0 (ms) 
  eGraded = -80 (mV)  
  vslopeGraded = 5.0 (mV)   
  vmidGraded = -40.0 (mV)   
  celltype = 1 
  alive = 1 
  cellid = -1 
  verbose = 0 
}


ASSIGNED {
  factor 
  eventflag 
  V (mV) 
  u (mV) 
  vpreGraded (mV)  
  gAMPA 
  gGraded  (uS) 
  I 
  delta 
  t0 
  sinfGraded  
}


INITIAL {
  V = vr
  u = 0.0
  t0 = t
  gAMPA = 0
  gGraded = 0
  I = 0
  delta = 0
  net_send(0,1) 
}



VERBATIM
char filename[1000]; 
ENDVERBATIM
PROCEDURE useverbose() { 
  VERBATIM
  #include<stdio.h>  
  verbose = (float) *getarg(1); 
  strcpy(filename, gargstr(2)); 
  ENDVERBATIM
}

STATE {
    sGraded
}


BREAKPOINT {
    SOLVE states METHOD cnexp
  delta = t-t0 

  
  gAMPA = gAMPA - delta*gAMPA/tauAMPA 
  gGraded = gbarGraded * sGraded
  
  factor = ((V+80)/60)*((V+80)/60)
  I = gAMPA*(V-0) + gGraded*(V - eGraded)
  
  
  Vpre = V
  V = V + delta*(k*(V-vr)*(V-vt) - u - I + Iin)/(C*100)  

  if (Vpre<=c && V>vpeak) {V=c+1} 

  
  if (celltype<5) {
    u = u + delta*a*(b*(V-vr)-u) 
  }
  else {
     
     if (celltype==5) {
       if (V<d) { 
        u = u + delta*a*(0-u)
       }
       else { 
        u = u + delta*a*((0.025*(V-d)*(V-d)*(V-d))-u)
       }
     }

     
     if (celltype==6) {
       if (V>-65) {b=0}
       else {b=15}
       u = u + delta*a*(b*(V-vr)-u) 
     }
     
     
     if (celltype==7) {
       if (V>-65) {b=2}
       else {b=10}
       u = u + delta*a*(b*(V-vr)-u) 
     }
  }

  t0=t 
  
  
  if (verbose>1) { 
    VERBATIM
    FILE *outfile; 
    outfile=fopen(filename,"a"); 
    fprintf(outfile,"%8.2f   cell=%6.0f   delta=%8.2f   gAMPA=%8.2f   gGraded=%8.2f   factor=%8.2f   I=%8.2f   V=%8.2f   u=%8.2f (timestep)\n",t,cellid,delta,gAMPA,gGraded,factor,I,V,u);
    fclose(outfile); 
    ENDVERBATIM
  }
}


NET_RECEIVE (wAMPA) {  
  INITIAL { wAMPA=wAMPA } 

  
  if (flag == 1) { 
    if (celltype < 4 || celltype == 5 || celltype == 7) { 
      WATCH (V>vpeak) 2 
    }
    else if (celltype == 4) { 
      WATCH (V>(vpeak-0.1*u)) 2 
    }
    else if (celltype == 6) { 
      WATCH (V>(vpeak+0.1*u)) 2 
    }
  } 
  
  
  else if (flag == 2) { 
    if (alive) {net_event(t)} 

    
    if (celltype < 4 || celltype == 7) {
      V = c 
      u = u+d 
    }
    
    else if (celltype == 4) {
      V = c+0.04*u 
      if ((u+d)<670) {u=u+d} 
      else {u=670} 
     }  
    
    else if (celltype == 5) {
      V = c 
     }  
    
    else if (celltype == 6) {
      V = c-0.1*u 
      u = u+d 
  }  
  
  gAMPA = 0 
  } 
  
  
  else {
    gAMPA = gAMPA + wAMPA
  }
  
  
  if (verbose>0) { 
    eventflag = flag
    VERBATIM
    FILE *outfile; 

    outfile=fopen(filename,"a"); 
    fprintf(outfile,"t=%8.2f   cell=%6.0f   flag=%1.0f   gAMPA=%8.2f    sinfGraded=%8.2f   gGraded=%8.2f   V=%8.2f   u=%8.2f (event)\n",t, cellid,eventflag,gAMPA,sinfGraded,gGraded,V,u);
    fclose(outfile); 

    ENDVERBATIM
}


}

DERIVATIVE states {
    sinfGraded = 1 / (1 + exp((vmidGraded - vpreGraded) / vslopeGraded))
    sGraded' = (sinfGraded - sGraded) / tauGraded
}