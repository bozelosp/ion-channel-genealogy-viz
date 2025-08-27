NEURON {
    POINT_PROCESS STDP 
    POINTER synweight 
    RANGE tauhebb, tauanti 
    RANGE hebbwt, antiwt 
    RANGE RLwindhebb, RLwindanti 
    RANGE useRLexp 
    RANGE RLlenhebb, RLlenanti 
    RANGE RLhebbwt, RLantiwt 
    RANGE wmax 
    RANGE softthresh 
    RANGE STDPon 
    RANGE RLon 
    RANGE verbose 
    RANGE tlastpre, tlastpost 
    RANGE tlasthebbelig, tlastantielig 
    RANGE interval 
    RANGE deltaw 
    RANGE newweight 
    RANGE skip 
}

ASSIGNED {
    synweight        
    tlastpre   (ms)    
    tlastpost  (ms)   
    tlasthebbelig   (ms)    
    tlastantielig  (ms)        
    interval    (ms)    
    deltaw
    newweight          
}

INITIAL {
    tlastpre = -1            
    tlastpost = -1           
    tlasthebbelig = -1      
    tlastantielig = -1  
    interval = 0
    deltaw = 0
    newweight = 0
}

PARAMETER {
    tauhebb  = 10  (ms)   
    tauanti  = 10  (ms)    
    hebbwt = 1.0
    antiwt = -1.0
    RLwindhebb = 10 (ms)
    RLwindanti = 10 (ms)
    useRLexp = 0   
    RLlenhebb = 100 (ms)
    RLlenanti = 100 (ms)
    RLhebbwt = 1.0
    RLantiwt = -1.0
    wmax  = 15.0
    softthresh = 0
    STDPon = 1
    RLon = 1
    verbose = 0
    skip = 0
}

NET_RECEIVE (w) {
    deltaw = 0.0 
    skip = 0
    
    if (verbose > 1)  { printf("t=%f (BEFORE) tlaspre=%f, tlastpost=%f, flag=%f, w=%f, deltaw=%f \n",t,tlastpre, tlastpost,flag,w,deltaw) }

    
    if ((flag == -1) && (tlastpre != t-1)) {   
        skip = 1 
        deltaw = hebbwt * exp(-interval / tauhebb) 
        if (softthresh == 1) { deltaw = softthreshold(deltaw) } 
        adjustweight(deltaw) 
        if (verbose > 1) { printf("Hebbian STDP event
        }

    
    else if ((flag == 1) && (tlastpost != t-1)) { 
        skip = 1 
        deltaw = antiwt * exp(interval / tauanti) 
        if (softthresh == 1) { deltaw = softthreshold(deltaw) } 
        adjustweight(deltaw) 
        if (verbose > 1) { printf("anti-Hebbian STDP event
        }


    
    if (skip == 0) {
        if (w >= 0) {           
            interval = tlastpost - t  
            if  ((tlastpost > -1) && (-interval > 1.0)) { 
                if (STDPon == 1) { 
                    if (verbose > 1) {printf("net_send(1,1)\n")}
                    net_send(1,1) 
                }
                if ((RLon == 1) && (-interval <= RLwindanti)) { tlastantielig = t } 
            }
            tlastpre = t 
        
        
        } else {            
            interval = t - tlastpre 
            if  ((tlastpre > -1) && (interval > 1.0)) { 
                if (STDPon == 1) { 
                    if (verbose > 1) {printf("net_send(1,-1)\n")}
                    net_send(1,-1) 
                }
                if ((RLon == 1) && (interval <= RLwindhebb)) { 
                    tlasthebbelig = t} 
            }
            tlastpost = t 
        }
    }
    if (verbose > 1)  { printf("t=%f (AFTER) tlaspre=%f, tlastpost=%f, flag=%f, w=%f, deltaw=%f \n",t,tlastpre, tlastpost,flag,w,deltaw) }
}

PROCEDURE reward_punish(reinf) {
    if (RLon == 1) { 
        deltaw = 0.0 
        deltaw = deltaw + reinf * hebbRL() 
        deltaw = deltaw + reinf * antiRL() 
        if (softthresh == 1) { deltaw = softthreshold(deltaw) }  
        adjustweight(deltaw) 
        if (verbose > 0) { printf("RL event
    }
}

FUNCTION hebbRL() {
    if ((RLon == 0) || (tlasthebbelig < 0.0)) { hebbRL = 0.0  } 
    else if (useRLexp == 0) { 
        if (t - tlasthebbelig <= RLlenhebb) { hebbRL = RLhebbwt } 
        else { hebbRL = 0.0 } 
    } 
    else { hebbRL = RLhebbwt * exp((tlasthebbelig - t) / RLlenhebb) } 
      
}

FUNCTION antiRL() {
    if ((RLon == 0) || (tlastantielig < 0.0)) { antiRL = 0.0 } 
    else if (useRLexp == 0) { 
        if (t - tlastantielig <= RLlenanti) { antiRL = RLantiwt } 
        else {antiRL = 0.0 } 
    }
    else { antiRL = RLantiwt * exp((tlastantielig - t) / RLlenanti) } 
}

FUNCTION softthreshold(rawwc) {
    if (rawwc >= 0) { softthreshold = rawwc * (1.0 - synweight / wmax) } 
    else { softthreshold = rawwc * synweight / wmax } 
}

PROCEDURE adjustweight(wc) {
   synweight = synweight + wc 
   if (synweight > wmax) { synweight = wmax }
   if (synweight < 0) { synweight = 0 }
}