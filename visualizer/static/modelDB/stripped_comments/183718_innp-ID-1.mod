NEURON {

    POINT_PROCESS InNp 
    NONSPECIFIC_CURRENT i
    RANGE delay, dur, per

    RANGE mean, stdev
    RANGE genevent, y1
    THREADSAFE 
    POINTER donotuse
}

UNITS {
    (nA) = (nanoamp)
    (mA) = (milliamp)
}

PARAMETER {
    delay (ms) 
    dur (ms) <0, 1e9> 
    per = 0.1 (ms) <1e-9, 1e9> 


    mean = 0 (nA)
    stdev = 1 (nA)
    genevent = 0 (1) 
      
      
}

ASSIGNED {
    on
    ival (nA)
    i (nA)
    donotuse
    t0 (ms)
    y0 (nA)
    y1 (nA)
}

INITIAL {
    on = 0
    ival = 0
    i = 0
    net_send(delay, 1)
}

PROCEDURE seed(x) {
    set_seed(x)
}

BEFORE BREAKPOINT {
    if (on==0) {
        i = 0
    } else {
        i = y0 + ((t-t0)/per)*(y1 - y0)
    }
}

BREAKPOINT {

}

FUNCTION yval() (nA) {

    yval = mean + nrand()*stdev 
}

NET_RECEIVE (w) {
    if (dur>0) {
        if (flag==1) {
            if (on==0) { 
                on=1
                net_send(dur,1) 
                net_send(per, 2) 
                t0 = t
                y0 = yval()
                y1 = yval()
                if (genevent==1) {
                    net_event(t) 
                }
            } else {
                if (on==1) { 
                    on=0
                    y0 = 0
                    y1 = 0
                }
            }
        }
        if (flag==2) {
            if (on==1) {
                net_send(per, 2) 
                t0 = t
                y0 = y1
                y1 = yval()
                if (genevent==1) {
                    net_event(t) 
                }
            }
        }
    }
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM





FUNCTION nrand() {
VERBATIM
    if (_p_donotuse) {
        


            _lnrand = nrn_random_pick(_p_donotuse);
    }else{
        
        if (_nt != nrn_threads) {
hoc_execerror("multithread random in InUnif"," only via hoc Random");
        }
ENDVERBATIM
        
        
        
        


        nrand = normrand(0,stdev/(1(nA)))
VERBATIM
    }
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
    void** pv = (void**)(&_p_donotuse);
    if (ifarg(1)) {
        *pv = nrn_random_arg(1);
    }else{
        *pv = (void*)0;
    }
 }
ENDVERBATIM
}