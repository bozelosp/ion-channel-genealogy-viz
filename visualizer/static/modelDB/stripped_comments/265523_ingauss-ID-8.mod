NEURON {
    POINT_PROCESS InGauss
    NONSPECIFIC_CURRENT i
    RANGE mean, stdev
    RANGE del, dur
    THREADSAFE 
    POINTER donotuse
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    del (ms) 
    dur (ms) <0, 1e9> 
    mean = 0 (nA)
    stdev = 1 (nA)
}

ASSIGNED {
    dt (ms)
    on
    per (ms)
    ival (nA)
    i (nA)
    donotuse
}

INITIAL {
    per = dt
    on = 0
    ival = 0
    i = 0
    net_send(del, 1)
}

PROCEDURE seed(x) {
    set_seed(x)
}

BEFORE BREAKPOINT {
    i = ival

}

BREAKPOINT { 
}

NET_RECEIVE (w) {
    if (dur>0) {
        if (flag==1) {
            if (on==0) { 
                on=1
                net_send(dur,1) 

                ival = stdev*grand() + mean 
                net_send(per, 2) 
            } else {
                if (on==1) { 
                    on=0
                    ival = 0
                }
            }
        }
        if (flag==2) {
            if (on==1) {
                ival = stdev*grand() + mean

                net_send(per, 2) 
            }
        }
    }
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM



FUNCTION grand() {
VERBATIM
    if (_p_donotuse) {
        


            _lgrand = nrn_random_pick(_p_donotuse);
    }else{
        
        if (_nt != nrn_threads) {
hoc_execerror("multithread random in InUnif"," only via hoc Random");
        }
ENDVERBATIM
        
        
        
        


        grand = normrand(0,1)

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