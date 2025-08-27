VERBATIM

#define SOBOL_MAXBIT 30
#define SOBOL_MAXDIM 6

float sobseq(int init)
{
  int j,k,l;
  unsigned long i,im,ipp;
  static float fac;
  static unsigned long in = 0,ix = 0, *iu[SOBOL_MAXBIT+1];
  static unsigned long mdeg[SOBOL_MAXDIM+1] = {0,1,2,3,3,4,4};
  static unsigned long ip[SOBOL_MAXDIM+1] = {0,0,1,1,2,1,4};
  static unsigned long iv[SOBOL_MAXDIM*SOBOL_MAXBIT+1] = {0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
  
  if (init) {
    if (iv[1] != 1)
      return -1.0;
    fac = 1.0/(1L << SOBOL_MAXBIT);
    for (j=1,k=0; j<=SOBOL_MAXBIT; j++,k+=SOBOL_MAXDIM)
      iu[j] = &iv[k];
    for (k=1; k<=SOBOL_MAXDIM; k++) {
      for (j=1; j<=mdeg[k]; j++)
	iu[j][k] <<= (SOBOL_MAXBIT-j);
      for (j=mdeg[k]+1; j<=SOBOL_MAXBIT; j++) {
	ipp = ip[k];
	i = iu[j-mdeg[k]][k];
	i ^= (i >> mdeg[k]);
	for (l=mdeg[k]-1; l>=1; l--) {
	  if (ipp & 1)
	    i ^= iu[j-l][k];
	  ipp >>= 1;
	}
	iu[j][k] = i;
      }
    }
    return 0;
  }
  im = in++;
  for (j=1; j<=SOBOL_MAXBIT; j++) {
    if (!(im & 1))
      break;
    im >>= 1;
  } 
  if (j > SOBOL_MAXBIT)
    fprintf(stderr, "SOBOL_MAXBIT (%d) too small in sobseq.\n", SOBOL_MAXBIT);
  im = (j-1)*SOBOL_MAXDIM;
  ix ^= iv[im+1];
  return ix*fac;
}

ENDVERBATIM

UNITS {
    (nA) = (nanoampere)
}

NEURON {
    POINT_PROCESS SobolPulses
    RANGE iPulse,iPI,F,delay,estimatedFreq, dur, amp, spkCount,gp,gi,maxPIcount
    ELECTRODE_CURRENT i
}

ASSIGNED {
    i                   (nA)
    tnext               (ms)
    tLastSpk            (ms)
    estimatedFreq       (1/s)
    count               (1)
    countPI             (1)
    erri                (1)
    errp                (1)
    cval                (nA)
    iPI                 (nA)
    iPulse              (nA)
}

PARAMETER {
    delay       = 1000  (ms)
    dur         = 0.5   (ms)
    amp         = 0.25  (nA)
    spkCount    = 10    (1)
    tau         = 0.01  (s)
    gp          = 0     (1)
    gi          = 0     (1)
    maxPIcount  = 3     (1)
    F           = 30    (1/s)
}   

INITIAL {
    tnext = -1
    VERBATIM
    sobseq(1);
    ENDVERBATIM
    erri = 0
    errp = 0
    cval = 0     (nA) 
    iPulse  = 0     (nA) 
    iPI     = 0     (nA) 
}

NET_RECEIVE(weight) { 
    LOCAL w,isi,frac
    INITIAL {
        count         = 0
        estimatedFreq = 0       
        tLastSpk      = 0
        countPI       = maxPIcount
    }
    if (tLastSpk > 0) {
        
        isi = (t - tLastSpk)/1000.0  
        w   = exp(-isi/tau)
        estimatedFreq = (1-w)/isi + w*estimatedFreq
        
        
        count = count + 1
        countPI = countPI + 1
    }
    if (t < tnext){
        tnext = -1
    }
    if (count >= spkCount && t > delay) {
        countPI = 0
        count   = 0
	frac    = sobseq(0)
        tnext   = t + frac*(1000.0/estimatedFreq)
        at_time(tnext)
	
    }
    if (countPI >= maxPIcount || count == (spkCount-1)){
        
        errp = F - estimatedFreq
        erri = erri + errp
        cval = (gp*errp + gi*erri)/1000.0
        
    }
    
	
    
    tLastSpk = t
}

BREAKPOINT {
    iPI = cval
    iPulse = 0
    if (t>tnext && t<=(tnext+dur)){
        iPulse = amp
        at_time(tnext+dur)
    }
    i=iPI+iPulse
    
}