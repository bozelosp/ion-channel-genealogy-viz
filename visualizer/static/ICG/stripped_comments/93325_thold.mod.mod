NEURON {
	POINT_PROCESS thold
	RANGE steadystate, reset, decaytc, lastspike, thresh, spikecount, prior2spike, burstmaxint, burstminsize, bursting, burstcount, burststart, burstno, freq
	POINTER nc_thresh
	}

UNITS {
	(mV) = (millivolt)
}

PARAMETER  {
	steadystate = -37.8 (mV)  
	reset = 200 (mV)	  
	decaytc = 11.5 (ms)	  
	nc_thresh (mV)		  
	lastspike (ms)		  
	thresh (mV)		  
	spikecount		  
	prior2spike		  
	burstmaxint = 1000
	burstminsize = 3
	bursting = 0
	burstcount = 0
	burststart = 0
	burstno = 0
	freq = 0
	}

ASSIGNED {
	v	(mV)
	}


BREAKPOINT {	
	thresh = steadystate + (reset - steadystate) * (exp((lastspike - t) / decaytc))
	nc_thresh = thresh  
}	

INITIAL {
	thresh = reset 
	nc_thresh = thresh
	lastspike = 0
	spikecount = 0
	prior2spike = -1
	bursting = 0
	burstcount = 0
	burststart = 0
	burstno = 0
	freq = 0
}


NET_RECEIVE(weight (microsiemens)) {
		if (spikecount > 0) {
		freq = (1/ (t - lastspike)* 1000) 
		if ((t - lastspike) < burstmaxint) {
			burstcount = burstcount + 1
			if (burstcount == 1) { burststart = lastspike }
			if (burstcount > burstminsize) { bursting = 1 }
		} 
		
		if ((t - lastspike) > burstmaxint) {
			bursting = 0
			burstcount = 0
			burststart = 0
		}
		}

		prior2spike = lastspike
		lastspike = t     
		spikecount = spikecount + 1  
		
	
}