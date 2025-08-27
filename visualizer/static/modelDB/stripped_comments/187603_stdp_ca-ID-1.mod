NEURON {
        POINT_PROCESS ExpSynCaSTDP
        RANGE tau, e, i, dd, dp, dtau, ptau, thresh, ca_thresh, wmax, wmin, D, cai
        NONSPECIFIC_CURRENT i
        USEION ca READ cai
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
}

PARAMETER {
        tau     = 3 (ms) <1e-9,1e9>
        e       = 0     (mV)
        dd      = 0.001 <0,1>   
        dp      = 0.00106       
        dtau    = 20 (ms)   
        ptau    = 20 (ms)   
        thresh  = -20 (mV)      
        ca_thresh  = 0.5 (mM)      
        wmax    = 0.001 (uS)
        wmin    = 0 (uS)
}

ASSIGNED {
        v (mV)
        cai (mM)
        i (nA)
    D
        tpost (ms)
}

STATE {
        g (uS)
}

INITIAL {
        g = 0
        D = 0
        tpost = -1e9
        net_send(0, 1)
}

BREAKPOINT {
        SOLVE state METHOD cnexp
        i = g*(v - e)
}

DERIVATIVE state {
        g' = -g/tau
}

NET_RECEIVE(w (uS), P, tpre (ms), wsyn) {
        INITIAL { P = 0 tpre = -1e9 wsyn = w}
        if (flag == 0) { 
	printf("entryca flag=%g t=%g w=%g A=%g tpre=%g tpost=%g g=%g ca=%g ca_thresh=%g\n", flag, t, w, P, tpre, tpost,g,cai,ca_thresh)
        printf("wsynca=%g\n", wsyn)
        
        g = g + wsyn
        P = P*exp((tpre-t)/ptau) + dp
        tpre = t
        wsyn = wsyn + wmax * D * exp((tpost-t)/dtau) 
        if (wsyn < wmin) { wsyn = wmin }
        }else if (flag == 2) { 
        printf("entryca flag=%g t=%g tpost=%g ca=%g\n", flag, t, tpost,cai)
        FOR_NETCONS(w1, P1, tpre1, wsyn1) {
            
            wsyn1 = wsyn1 + wmax * P1 * exp(-(t - tpre1)/ptau) 
            if (wsyn1 > wmax) { wsyn1 = wmax }
        }
        D = D*exp((tpost-t)/dtau) - dd
        
        
        
        tpost = t
        } else { 
        printf("entryca flag=%g t=%g ca=%g\n", flag, t,cai)
                WATCH (cai > ca_thresh) 2
        }
}