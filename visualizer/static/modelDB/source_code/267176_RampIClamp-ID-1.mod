    COMMENT
    Point process for generation of ascending and descending current injection over time.
    ENDCOMMENT
  
    NEURON {
            POINT_PROCESS RampIClamp
            RANGE  dur, amp,i   
            ELECTRODE_CURRENT i
    }

    UNITS {
            (nA) = (nanoamp)
          }

    ASSIGNED {
        i(nA)   
        tc2 (ms)
        tc3 (ms)
	dur[3] (ms)
	amp[3] (mV)
    }

     INITIAL {
       tc2 = dur[0] + dur[1]
        tc3 = tc2 + dur[2]
       
}


    BREAKPOINT {
          
    if (t < dur[0]) {
                i = amp[0]
        }else if (t < tc2) {
                i =((amp[1]-amp[0])/dur[1])*(t-dur[0])+amp[0]  
        }else if (t < tc3) {
                i = ((amp[2]-amp[1])/dur[2]) *(t-tc2)+amp[1]   
        }else {
                i= amp[2]
                 
        }
      }