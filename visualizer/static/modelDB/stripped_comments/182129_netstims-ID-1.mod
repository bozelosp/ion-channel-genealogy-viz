NEURON	{ 
  POINT_PROCESS NetStims
  RANGE y
  RANGE interval, number, start
  RANGE noise
  RANGE freqhz
  RANGE q,prob
}

PARAMETER {
	interval	= 10 (ms) <1e-9,1e9>
	number	= 1000000 <0,1e9>	
	start		= 50 (ms)	
	noise		= 0 <0,1>	
	freqhz=10
	q=0.30
	prob=0.50
	
}

ASSIGNED {
	y
	event (ms)
	on
	end (ms)
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	y = 0
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
	
		event = start
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = t
		end = t + 1e-6 + invl(interval)*(number-1)
	}
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 
	}
	if (noise == 0) {
		invl = mean
	}else{
        invl = twoRndGamma(freqhz,q, prob)
		
		
		
		
	}
}

PROCEDURE event_time() {
	if (number > 0) {
		event = event + invl(interval)
	}
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
		y = 2
		net_event(t)
		event_time()
		if (on == 1) {
			net_send(event - t, 1)
		}
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = 0
	}
}

FUNCTION twoRndGamma(freqhz,q, prob) {
    LOCAL rnd,time,tt,ts,tl
	
	tt=1000/freqhz
	
	ts=tt/prob
	
	if (q==1){
	tl=tt-q*ts
	}else{
	tl=(tt-q*ts)/(1-q)
	
	
	}
	
	rnd=scop_random()
	
	if (rnd<q) {
	    time=ts
	  }else{
	  time=tl
	  }
    rnd=scop_random()
	if (q==0) {twoRndGamma=time
	} else{	
	twoRndGamma = -time*log(1-rnd)}

}