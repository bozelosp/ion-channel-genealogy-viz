INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON	{ 
  POINT_PROCESS gen
  GLOBAL seed
  RANGE x, spk, lastspk
  RANGE fast_invl, slow_invl, burst_len, start, end
  RANGE noise
}

PARAMETER {
	fast_invl	= 1		
	slow_invl	= 50		
	burst_len	= 10		
	start		= 50		
	end		= 1e10		
	noise		= 0		
	seed            = 53892         
}

ASSIGNED {
	scntr
	lcntr
	bcntr
	burst
	x
        spk
        lastspk
	dt
}

INCLUDE "presyn.inc"

INITIAL {
	burst   = 0
	scntr	= fast_invl
	bcntr	= burst_len
	if (noise != 0) {
	  set_seed(seed)
	  scntr	= noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
	  bcntr	= noise * fpoisrand(burst_len) + (1 - noise) * burst_len
	}
	lcntr	= start - fast_invl
	x = -90
	spk = -90
}	

BREAKPOINT {
  SOLVE generate
}

PROCEDURE generate() {	

  if (t < end) {

    
    if (t > lastspk+1) { spk = -90 }
    x = -90
    scntr = scntr - dt
    lcntr = lcntr - dt	

    if (burst) {
      

      if (scntr <= dt) {
	
	
	if (bcntr <= 1) {	
	  

	  burst = 0
	  if (noise==0) {
	    lcntr = slow_invl
	  } else {
	    lcntr = noise * fpoisrand(slow_invl) + (1 - noise) * slow_invl
	  }
	}

	x = 50		
        spk = 50
        lastspk = t
	newspike()

	bcntr = bcntr - 1
	if (noise==0) {
	  scntr = fast_invl
	} else {
	  scntr = noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
	}
      } 	

    } else {
      

      if (lcntr <= dt) {
	
	if (noise==0) {
	  bcntr = burst_len
	} else {
	  bcntr = noise * fpoisrand(burst_len) + (1 - noise) * burst_len
	}
	burst = 1
	
      }

    }	

    VERBATIM
    return 0;
    ENDVERBATIM
  }
}	

FUNCTION fgauss(x,mean,std_dev) {
	fgauss = gauss(x,mean,std_dev)
}

FUNCTION fpoisrand(mean) {
  if (mean > 700) {
    fpoisrand = 4. * poisrand(mean/4.)  
  } else {
    fpoisrand = poisrand(mean)
  }
}