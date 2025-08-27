NEURON {
            POINT_PROCESS SinClamp
            RANGE delay, dur, pkamp, freq, phase, bias
            ELECTRODE_CURRENT i
    }

    UNITS {
            (nA) = (nanoamp)
                 }

    PARAMETER {
            delay=5   (ms)
            dur=200   (ms)
            pkamp=0 (nA)
            freq=1  (Hz)
            phase=0 (rad)
            bias=0  (nA)
            PI=3.14159265358979323846
    }

    ASSIGNED {
            i (nA)
    }

    BREAKPOINT {
           at_time(delay)
           at_time(delay + dur)

           if (t < delay) {
             i=0   
       	   } else {
                if (t < delay+dur) {
	        	i = pkamp*sin(2*PI*freq*(t-delay)/1000+phase)+bias
		        
         	    } else { 
	               i = 0
  		    }
	          }
    	   }