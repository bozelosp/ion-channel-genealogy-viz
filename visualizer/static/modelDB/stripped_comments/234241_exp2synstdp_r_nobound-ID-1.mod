NEURON {
        POINT_PROCESS Exp2Syn_Inh_STDP_nobound
        RANGE tau1, tau2, e, i, dd, dp, dtau, ptau, thresh, wmax, wmin, D, shift, g, mean, std, learning_rate
        NONSPECIFIC_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
}

PARAMETER {
        tau1     = 0.1 (ms) <1e-9,1e9>
        tau2    = 10 (ms) <1e-9,1e9>
        e       = -73     (mV)
        dd      = 0.001 <0,1>   
        dp      = 0.001       
        dtau    = 20 (ms)   
        ptau    = 20 (ms)   
        thresh  = -20 (mV)      
        wmax    = 0.001 (uS)
        wmin    = 0 (uS)
        shift   = 0 (ms) 	
        mean    = 0 (ms)
        std     = 0 (ms)
        learning_rate = 1

}

ASSIGNED {
        v (mV)
        i (nA)
        D
        tpost (ms)
        g        
        factor
        jitter
}

STATE {
        A (uS)
        B (uS)
}

INITIAL {
        LOCAL tp
        A = 0
        B = 0
        D = 0
        tpost = -1e9
        jitter = normrand(mean, std)
        net_send(0, 1)
        if (tau1/tau2 > .9999){
            tau1 = .9999*tau2
        }
        tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
        factor = -exp(-tp/tau1) + exp(-tp/tau2)
        factor = 1/factor
}

BREAKPOINT {
        SOLVE state METHOD cnexp
        g = B-A
        i = g*(v - e)
}

DERIVATIVE state {
        A' = -A/tau1
        B' = -B/tau2
}

NET_RECEIVE(w (uS), P, tpre (ms), wsyn) {
        INITIAL { P = 0 tpre = -1e9 wsyn = w}

        if (flag == 0) { 
                
                
                
        net_send(jitter,5)
        net_send(jitter+shift, 3)        
        
        
        
        
        
        }else if (flag == 2) { 
                
        FOR_NETCONS(w1, P1, tpre1, wsyn1) {
            wsyn1 = wsyn1 + learning_rate * wmax * P1 * exp(-(t-tpre1)/ptau) 
            if (wsyn1 < wmin) { wsyn1 = wmin }
        }
        
        D = D*exp((tpost-t)/dtau) + dd
        tpost = t

        } else if (flag == 3) { 
                
        P = P*exp((tpre-t)/ptau) - dp
        tpre = t
        wsyn = wsyn + learning_rate * wmax * D * exp((tpost-t)/dtau) 
        

        
        
        
        
        
        
        
        
        } else if (flag == 5) {
        A = A + wsyn*factor
        B = B + wsyn*factor
        } else { 
        
                WATCH (v > thresh) 2
        }
}