NEURON {
  POINT_PROCESS Izhi2007
  RANGE C, k, vr, vt, vpeak, a, b, c, d, Iin, Vpre, tauAMPA, tauNMDA, tauGABAA, tauGABAB, tauOpsin, celltype, alive, cellid, verbose, t0
}


UNITS {
  (mV) = (millivolt)
  (uM) = (micrometer)
}


PARAMETER {
  C = 100 
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
  tauAMPA = 5 (ms) 
  tauNMDA = 150 (ms) 
  tauGABAA = 6 (ms) 
  tauGABAB = 150 (ms) 
  tauOpsin = 50 (ms) 
  celltype = 1 
  alive = 1 
  cellid = -1 
  verbose = 0 
}



ASSIGNED {
  
  factor 
  eventflag 

  
  

  V 
  u 
  gAMPA 
  gNMDA 
  gGABAA 
  gGABAB 
  gOpsin 
  I 
  delta 
  t0 
}



INITIAL {
  V = vr
  u = 0.2*V
  t0 = t
  gAMPA = 0
  gNMDA = 0
  gGABAA = 0
  gGABAB = 0
  gOpsin = 0
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



BREAKPOINT {
  delta = t-t0 

  
  gAMPA = gAMPA - delta*gAMPA/tauAMPA 
  gNMDA = gNMDA - delta*gNMDA/tauNMDA 
  gGABAA = gGABAA - delta*gGABAA/tauGABAA 
  gGABAB = gGABAB - delta*gGABAB/tauGABAB 
  gOpsin = gOpsin - delta*gOpsin/tauOpsin 
  
  
  factor = ((V+80)/60)*((V+80)/60)
  I = gAMPA*(V-0) + gNMDA*factor/(1+factor)*(V-0) + gGABAA*(V+70) + gGABAB*(V+90) + gOpsin*(V-0) 
  
  
  Vpre = V
  V = V + delta*(k*(V-vr)*(V-vt) - u - I + Iin)/C  

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
        u = u + delta*a*((0.025*(V-d)^3)-u)
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
    fprintf(outfile,"%8.2f   cell=%6.0f   delta=%8.2f   gAMPA=%8.2f   gNMDA=%8.2f   gGABAA=%8.2f   gGABAB=%8.2f   gOpsin=%8.2f   factor=%8.2f   I=%8.2f   V=%8.2f   u=%8.2f (timestep)\n",t,cellid,delta,gAMPA,gNMDA,gGABAA,gGABAB,gOpsin,factor,I,V,u);
    fclose(outfile); 
    ENDVERBATIM
  }
}


NET_RECEIVE (wAMPA, wNMDA, wGABAA, wGABAB, wOpsin) {  
  INITIAL { wAMPA=wAMPA wNMDA=wNMDA wGABAA=wGABAA wGABAB=wGABAB wOpsin=wOpsin} 

  
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
    gNMDA = 0
    gGABAA = 0
    gGABAB = 0
    gOpsin = 0
  } 
  
  
  else {
    gAMPA = gAMPA + wAMPA
    gNMDA = gNMDA + wNMDA
    gGABAA = gGABAA + wGABAA
    gGABAB = gGABAB + wGABAB
    gOpsin = gOpsin + wOpsin
  }
  
  
  if (verbose>0) { 
    eventflag = flag
    VERBATIM
    FILE *outfile; 

    outfile=fopen(filename,"a"); 
    fprintf(outfile,"t=%8.2f   cell=%6.0f   flag=%1.0f   gAMPA=%8.2f   gNMDA=%8.2f   gGABAA=%8.2f   gGABAB=%8.2f   gOpsin=%8.2f   V=%8.2f   u=%8.2f (event)\n",t, cellid,eventflag,gAMPA,gNMDA,gGABAA,gGABAB,gOpsin,V,u);
    fclose(outfile); 

    ENDVERBATIM
  }
  
  
}