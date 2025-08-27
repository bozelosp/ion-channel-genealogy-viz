    COMMENT
    Point process for generation of ascending and descending current injection over time.

    Three different stimulation profiles were used in the network model validation with
        the following parameters:
        ALL - Random noise with min = 6 (nA) and std = 0.2 (nA)
        (1) Ramp-up and hold profile
            dur = 3000 (ms)
            ramp_dur = 1000 (ms)
        (2) Ramp-up and hold profile with longer ramp
            dur = 4000 (ms)
            ramp_dur = 2000 (ms)
        (3) Ramp-up and ramp-down profile
            dur = 4000 (ms)
            ramp_dur = 2000 (ms)
    
    All three stimulation profiles were ran at three different peak amplitudes
        (10.5, 14.8, and 20 nA) corresponding to 10%, 40%, and 75% MVC. 
    ENDCOMMENT

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
                :: For ramp-up and hold stimulation profile
                amp = pkamp

                :: For ramp-up and ramp-down stimulation profile
                ::  (uncomment below and comment previous profile in else statement)
                :amp = -((pkamp-min)*((t-dur)/ramp_dur)) + min
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
                : turn it on
                on = 1
                : prepare to turn it off
                net_send(dur,1)
            } else {
                : turn it off
                on = 0
            }
        }
    }