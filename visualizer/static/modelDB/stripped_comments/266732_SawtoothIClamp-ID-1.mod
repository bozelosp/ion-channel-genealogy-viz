NEURON {
            POINT_PROCESS SawtoothIClamp
            RANGE del, tp, pkamp, bias
            ELECTRODE_CURRENT i
    }

    UNITS {
            (nA) = (nanoamp)
                 }

    PARAMETER {
            del=0   (ms)
            tp=10000   (ms)
            pkamp=10 (nA)
			bias=0 (nA)
			pie=3.14159
    }

    ASSIGNED {
            i (nA)
    }

    BREAKPOINT {
           at_time(del)

           if (t < del) {
          i=0   
       }else{ 
          i= pkamp/pie*acos(sin(2*pie/tp*(t+tp/4-del)))+bias  
	}}