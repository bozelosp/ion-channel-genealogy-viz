NEURON {
            POINT_PROCESS RampIClampNoise
            RANGE del, dur, pkamp, bias, ramp_dur,std
            ELECTRODE_CURRENT i
    }

    UNITS {
            (nA) = (nanoamp)
          }

    PARAMETER {
            del=0   (ms)
            dur=3000   (ms)
            ramp_dur = 1000 (ms)
            pkamp=20 (nA)
            bias=0 (nA)
			min = 6 (nA)
            std=0.2 (nA)
    }

    ASSIGNED {
            i (nA)
            ival (nA)
            amp (nA)
            noise (nA)
            on (1)
    }

    INITIAL {
        i = 0
        on = 0
        net_send(del,1)
    }

    PROCEDURE seed(x) {
        set_seed(x)
    }

    BEFORE BREAKPOINT {
        if (on) {
            noise = normrand(0,std*1(/nA)*1(nA))
            if (t < ramp_dur) {
                amp = ((pkamp-min)*(t/ramp_dur)) + min
            } else{
                
                amp = pkamp

                
                
                
            }
            ival = amp + noise
        } else {
            ival = 0
        }
    }

    BREAKPOINT {
        i = ival
    }

    NET_RECEIVE (w) {
        if (flag == 1) {
            if (on == 0) {
                
                on = 1
                
                net_send(dur,1)
            } else {
                
                on = 0
            }
        }
    }