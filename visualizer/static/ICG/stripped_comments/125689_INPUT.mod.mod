UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    NONSPECIFIC_CURRENT i
    SUFFIX INPUT

    GLOBAL  spikedur
        RANGE   Spike, lastspk, spk

    GLOBAL seed
    GLOBAL spkgenmode
        RANGE fast_invl, slow_invl, burst_len, start, end
        RANGE noise
        RANGE on_times

}


PARAMETER {

    v
    i

    spikedur = 1.0  (ms)
        refact   = 1.5  (ms)

        spkgenmode      = 0             
        fast_invl       = 1             
        slow_invl       = 100           
        burst_len       = 1             
        start           = 50            
        end             = 200           
        noise           = 0             
        seed            = 53892         

}

ASSIGNED {

        Spike
        spk
        lastspike

        pulsecntr
        scntr
        lcntr
        bcntr
        burst
        dt
        on_times[10]            
                                                

}

INITIAL {
        Spike = 0
        lastspike = -9e4

        pulsecntr = -1
        burst   = 0
        scntr   = fast_invl
        bcntr   = burst_len
        if (noise != 0) {
          set_seed(seed)
          scntr = noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
          bcntr = noise * fpoisrand(burst_len) + (1 - noise) * burst_len
        }

        lcntr   = start
}

BREAKPOINT {
        SOLVE generate
        i = Spike*(v+20)     
}



PROCEDURE generate() {
  if (t < end) {
    
    if (t > lastspike+1) { spk = 0 }
    Spike = 0
    scntr = scntr - dt
    lcntr = lcntr - dt
    if (burst) {
      
      if (scntr <= dt+dt/2) {
                        
                        if (bcntr <= 1) {
                                
                                burst = 0
                                if (spkgenmode==0) {
                                if (noise==0) {
                                lcntr = slow_invl
                            } else {
                            lcntr = noise * fpoisrand(slow_invl) + (1 - noise) * slow_invl
                                }
                                } else if  (spkgenmode==1) {
                                lcntr = on_times[pulsecntr] + dt
                                }
                }

                        Spike = 1
                        spk = 1
         lastspike = t

                        bcntr = bcntr - 1
                        if (noise==0) {
                                scntr = fast_invl + dt
                        } else {
                                scntr = noise * fpoisrand(fast_invl) + (1 - noise) * fast_invl
                        }
        }                                                                                       
        } else {                                                                        
      
      if (lcntr <= dt+dt/2) {                           
                        
                if (noise==0) {
                        bcntr = burst_len
                } else {
                        bcntr = noise * fpoisrand(burst_len) + (1 - noise) * burst_len
                }
                burst = 1
        pulsecntr = pulsecntr + 1
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