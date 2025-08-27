NEURON {
 POINT_PROCESS SpikeStats
 RANGE binsnum, ignoreBefore, stimPeriod, poissDt, sinsum, cossum
 POINTER n, vs, phase, pst, sinphs, cosphs
}

UNITS {
 PI = (pi) (1)
}

CONSTANT {radToDeg = 57.295779513082322865} 

PARAMETER {
 binsnum = 80                 <0,1e9>
 ignoreBefore = 0 (ms)        <0,1e9>
 stimPeriod  = 1 (ms)         <0,1e9>
 poissDt = 0.0125 (ms)        <1e-9,1e9>
 
 
 
 
}

ASSIGNED {
 sinsum
 cossum
 n      
 vs     
 phase  
 pst    
 sinphs 
 cosphs 
}

INITIAL {
 LOCAL phsi
 VERBATIM
  
#define PST    (&pst)
#define SINPHS (&sinphs)
#define COSPHS (&cosphs)
#define I_INT  ((int)_li)
#define PHSI   (_lphsi)
 ENDVERBATIM
 sinsum = 0
 cossum = 0
 if (pointersValid()) {
  n = 0
  vs = 0
  phase = -1e-9
  FROM i=0 TO binsnum-1 {
   phsi = 2*PI/binsnum*i
   VERBATIM
    PST[I_INT]    = 0.0;       
    SINPHS[I_INT] = sin(PHSI); 
    COSPHS[I_INT] = cos(PHSI); 
   ENDVERBATIM
  }
 }
}

NET_RECEIVE(change) {
 LOCAL i
 VERBATIM
  
#define PST    (&pst)
#define SINPHS (&sinphs)
#define COSPHS (&cosphs)
#define I_INT  ((int)_li)
 ENDVERBATIM
 if (t > ignoreBefore) {
  n = n + 1
  UNITSOFF 
   i = fmod(floor(fmod(t,stimPeriod)/poissDt),binsnum)
  UNITSON
  VERBATIM
   (PST[I_INT])++ ;         
   sinsum += SINPHS[I_INT]; 
   cossum += COSPHS[I_INT]; 
  ENDVERBATIM
  vs = sqrt(sinsum*sinsum + cossum*cossum)/n
  phase = atan2(sinsum,cossum)*radToDeg
  if (phase < 0) { phase = phase + 360 }
 }
}

FUNCTION pointersValid() {
 if ((!nrn_pointing(n))||(!nrn_pointing(vs))||(!nrn_pointing(phase))||
     (!nrn_pointing(pst))||(!nrn_pointing(sinphs))||
     (!nrn_pointing(cosphs))) {
  printf("SpikeStats
  printf("must all be set using setpointer!\n")
  pointersValid = 0
 } else {
  pointersValid = 1
 }
}