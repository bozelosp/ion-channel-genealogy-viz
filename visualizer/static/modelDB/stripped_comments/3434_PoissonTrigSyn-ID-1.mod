NEURON {
 POINT_PROCESS PoissonTrigSyn
 RANGE tau, refraction, gmax, e, i, g, k, firing
 POINTER prob, randseed
 NONSPECIFIC_CURRENT i
}

UNITS {
 (nA) = (nanoamp)
 (mV) = (millivolt)
 (uS) = (microsiemens)
}

PARAMETER {
 tau=.1 (ms)
 refraction = 1 (ms)
 gmax=0  (uS)
 e=0 (mV)
 randseed
}

ASSIGNED {
 i (nA)
 g (uS)
 k (/ms)
 firing
 prob
 v (mV)
}

CONSTANT { exp1 = 2.7182818284590452354} 

STATE { A (uS) G (uS) }

INITIAL {
 k = 1/tau
 A = 0
 G = 0
 firing = 0
 if (!nrn_pointing(prob)) {
  printf("PoissonTrigSyn
 }
 if (!nrn_pointing(randseed)) {
  printf("PoissonTrigSyn
 }
}


BREAKPOINT {
 SOLVE state METHOD sparse
 i = G*(v - e)
 g = G
}

KINETIC state {
 ~ A <-> G (k, 0)
 ~ G ->    (k)
}

NET_RECEIVE (weight) {
 if (flag == 0) {  
  if ((prob > unitRand()) && (!firing)){
   state_discontinuity(A, A + gmax*exp1)
   firing = 1
   if (refraction >= 0) {
    net_send(refraction, firing) 
   } else {
    firing = 0                   
   }
  }
 } else {
  firing = 0  
 }
}

FUNCTION getseedCD() {
 VERBATIM
#define RANDSEED_UINT (*((unsigned int*)(&randseed)))
#define GETSEEDCD (_lgetseedCD)
  GETSEEDCD = (double)RANDSEED_UINT; 
 ENDVERBATIM
}

PROCEDURE initseedCD(seed) {
 VERBATIM
#define RANDSEED_UINT (*((unsigned int*)(&randseed)))
#define SEED (_lseed)
  RANDSEED_UINT = (int)SEED;   
   
 ENDVERBATIM
}

FUNCTION unitRand() { 
 VERBATIM
#define RANDSEED_UINT (*((unsigned int*)(&randseed)))
#define UNITRAND (_lunitRand)
   
   RANDSEED_UINT = RANDSEED_UINT * 1103515245 + 12345;
   
   UNITRAND = ((RANDSEED_UINT >> 16) & 0x7FFF)/(((unsigned int)(0x7FFF))+1.0);
 ENDVERBATIM
}