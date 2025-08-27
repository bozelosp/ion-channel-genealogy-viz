NEURON {
    POINT_PROCESS synph
    RANGE nsyn, e1del, e1flag, e2del, e2flag
    RANGE ratio, noampablock, nonmdablock, sloc
}

PARAMETER {
    
    nsyn = 1

    
    
    e1del = 2 (ms)
    e1flag = 1
    e2del = 22 (ms)
    e2flag = 0

    
    ratio = 1           
    noampablock = 1
    nonmdablock = 1

    
    
    sloc = -1
}