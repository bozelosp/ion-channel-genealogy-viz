NEURON  {
    ARTIFICIAL_CELL GammaStim
    RANGE interval, start, duration, order, noise, refractoryPeriod
}

PARAMETER {
    interval = 10 (ms) <1e-9,1e9>	
    start = 1 (ms)       		    
    noise = 0 <0,1>       		    
    duration = 1000 (ms)		    
    order = 1 <1,6>                 
    refractoryPeriod = 0 (ms)
}

ASSIGNED {
    event (ms)
    on
    end (ms)
}

PROCEDURE seed(x) {
    set_seed(x) 
                
}

INITIAL {

    on = 0 
    if (order < 1 || order > 6) {
        order = 1
    }
    if (noise < 0) {
        noise = 0
    }
    if (noise > 1) {
        noise = 1
    }
    if (start >= 0) {
        
        
        event = start + invl(interval) - interval*(1. - noise)
        
        if (event < 0) {
            event = 0
        }
        net_send(event, 3)
    }
}

PROCEDURE init_sequence(t(ms)) {
    on = 1
    event = t
    end = t + 1e-6 + duration
}

FUNCTION invl(mean (ms)) (ms) {

    

    if (mean <= 0.) {
        mean = .01 (ms)
    }
    if (noise == 0) {
        invl = mean
    }else{
        invl = (1. - noise)*mean + noise*meanRndGamma(order, refractoryPeriod, mean)
    }
}

PROCEDURE event_time() {
    event = event + invl(interval)
    if (event > end) {
        on = 0
    }
}

NET_RECEIVE (w) {
    if (flag == 0) { 
        if (w > 0 && on == 0) { 
            init_sequence(t)
            net_send(0, 1)
                    
        }else if (w < 0 && on == 1) { 
            on = 0
        }
    }
    if (flag == 3) { 
        if (on == 0) {
            init_sequence(t)
            net_send(0, 1)
        }
    }
    if (flag == 1 && on == 1) {
        net_event(t) 
        event_time()
        if (on == 1) {
            net_send(event - t, 1)
        }
        net_send(.1, 2)
    }
}

FUNCTION meanRndGamma(gammaOrder(1), refractoryPeriod(ms), mean(ms)) (1) {

    

	LOCAL x

	x = 1.0
	FROM i = 0 TO gammaOrder-1 {
	    x = x * scop_random()
    }
	x = -log(x) * (interval - refractoryPeriod) / gammaOrder
	meanRndGamma = x + refractoryPeriod
}