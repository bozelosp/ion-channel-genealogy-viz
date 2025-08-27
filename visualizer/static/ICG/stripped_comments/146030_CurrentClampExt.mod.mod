NEURON {

    POINT_PROCESS CurrentClampExt
    RANGE del, dur, amp, repeat, i
    ELECTRODE_CURRENT i
    
}

UNITS {

    (nA) = (nanoamp)
    
}


PARAMETER {

    del (ms)
    dur (ms)    <0,1e9>
    amp (nA)
    repeat = 0
    
}

ASSIGNED { 
    
    i (nA) 
    
}

INITIAL {

    i = 0
    
}

BREAKPOINT {

    LOCAL shiftedTime, beginNextCycle
    
    shiftedTime = t
    
    at_time(del)
    at_time(del+dur)
    
    if (repeat == 1)
    {
        beginNextCycle = 0
        
        while (shiftedTime > del + dur)
        {
            shiftedTime = shiftedTime - (del + dur)
            beginNextCycle = beginNextCycle + (del + dur)
        }
        
        at_time(beginNextCycle + del)       ? to inform CVODE about a future discontinuity
        at_time(beginNextCycle + del+dur)
    }

    if (shiftedTime < del + dur && shiftedTime >= del) {
        i = amp
    }else{
        i = 0
    }
    
    
}