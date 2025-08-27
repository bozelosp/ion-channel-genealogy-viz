NEURON {
	POINT_PROCESS noisesyn
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	POINTER ptr
	RANGE g,start,spikedur,spikefreq,weight,nospike_tau,spike_tau,normalmean,normalstd,poisson_mean,gid,syn_index
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.05 (ms) <1e-9,1e9>
	tau2 = 5.3 (ms) <1e-9,1e9>
	e=0	(mV)
	
	start = 10 (ms)
	nospike_tau=3.3333 (ms) 
	spike_tau = 0.6667 (ms) 
	spikedur = 150 (ms)
	spikefreq = 2 (hz)
	normalmean = 0 (ms)
	normalstd = 4.4721 (ms) 
	weight = 0.00053407075(uS) 
	poisson_mean = 0.8 
	gid = 0 
	syn_index = 0 
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	t_master (ms) 
	interspike_int (ms)
	t_out (ms) 
	temp_time (ms) 
	ptr
	cachedNormal 
	haveCachedNormal 
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1) 
	factor = -exp(-tp/tau1) + exp(-tp/tau2) 
	factor = 1/factor 
	
	setrand(gid,syn_index) 
	cachedNormal=0 
	haveCachedNormal=0 
	
	t_master=start
	interspike_int = 1000/spikefreq-spikedur 
	net_send(t_master + normal(normalmean,normalstd),1) 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(dummy (uS)) {
	if(flag==1) { 
		chweight(my_poisson(poisson_mean)*weight*factor) 							
		t_master=t_master+spikedur 							
		t_out=t_master + normal(normalmean,normalstd) 		
		
		net_send(negexp(spike_tau),3) 						
	} else if(flag==2) { 
		
		t_master=t_master+interspike_int 					
		t_out=t_master + normal(normalmean,normalstd)
		temp_time=negexp(nospike_tau)  
		if(t + temp_time<t_out) { 
			net_send(temp_time,4) 					
		} else{ 
			net_send(t_out-t,1) 
		}
	} else if(flag==3) { 
		
		if(t >= t_out) { 
			net_send(0,2) 
		} else{ 
			chweight(my_poisson(poisson_mean)*weight*factor) 
			net_send(negexp(spike_tau),3)					
		}	
	} else if(flag==4) { 
		
		
		
		
		
		chweight(my_poisson(poisson_mean)*weight*factor) 						
		temp_time=negexp(nospike_tau)  
		if(t + temp_time<t_out) { 
			net_send(temp_time,4) 
		} else{ 
			net_send(t_out-t,1)					
		}
	}
}



VERBATIM
#define VOIDCAST void** vp = (void**)(&(_p_ptr))
extern void * nrnran123_newstream(int,int);
extern void nrnran123_deletestream(void *);
extern double nrnran123_dblpick(void *);
ENDVERBATIM

PROCEDURE setrand(id1,id2) {
	VERBATIM
	VOIDCAST;
	if(*vp) {
		nrnran123_deletestream(*vp);
	} 
	*vp = nrnran123_newstream((int) _lid1,(int) _lid2);
	ENDVERBATIM
} 

FUNCTION pick() {
	VERBATIM
	VOIDCAST;
	_lpick = nrnran123_dblpick(*vp);
	ENDVERBATIM
}

PROCEDURE chweight(delta) {
	
	A = A + delta
	B = B + delta
}

FUNCTION normal(pMean,pStdDev) { 
VERBATIM
if (haveCachedNormal == 1) {
	haveCachedNormal = 0;
	_lnormal = cachedNormal * _lpStdDev + _lpMean ;
} else {

		for(;;) {
			double u1 = pick();
			double u2 = pick();
			double v1 = 2 * u1 - 1;
			double v2 = 2 * u2 - 1;
			double w = (v1 * v1) + (v2 * v2);
 	    
	
	
	
	
			if (w <= 1) {
				double y = sqrt( (-2 * log(w)) / w);
				double x1 = v1 * y;
				double x2 = v2 * y;
	
				haveCachedNormal = 1;
				cachedNormal = x2;
				_lnormal = x1 * _lpStdDev + _lpMean;
				return _lnormal;
			}
		}
    }
ENDVERBATIM
}

FUNCTION negexp(mean) { 
	negexp = -mean*log(pick())
}

FUNCTION my_poisson(mean) { 
	LOCAL bound,count,product
	bound = exp(-1.0 * mean)
	count = 0
	product = 1
	while(product >= bound) {
		product = product * pick()
		count = count + 1
	}
	my_poisson = count - 1
}