NEURON	{ 
  POINT_PROCESS GaussStim
  RANGE x
  RANGE interval, number, start,factor,refrac
  
  RANGE N_backward,N_forward,N_normal,N_total
  RANGE rnd
  RANGE rand
  
}

PARAMETER {
	interval	= 10 (ms) <1e-9,1e9>
	number	= 10000 <0,1e9>	
	start		= 10 (ms)	
	factor =4  
	refrac		=0.5 <0,3>    
							
}

ASSIGNED {
	x
	event (ms)	
	on
	end (ms)
	m (ms)            
	diff (ms)
	N_forward  
	N_backward  
	N_normal
	N_total
	rnd
	rand
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	x = 0
	diff=0
	m=interval/2  
	N_forward=0
	N_backward=0
	N_normal=0
	N_total=0
	if (start >= 0 && number > 0) {
		
		
		event = start + invl(m) 
		
		
		if (event < 0) {
			event = 0
		}
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = t
		end = t + 1e-6 + interval*(number-1)
		
	}
}

FUNCTION invl(mean (ms)) (ms) { LOCAL std
	if (mean <= 0.) {
		mean = .01 (ms) 
	}
	std=mean/factor  
	invl = normrand(mean, std)  
	
	if(invl>=interval) { 
		invl=fmod(invl,interval)
		N_forward=N_forward+1
		}else if(invl<0) { 

			invl=fmod(invl,interval)+interval
			
			N_backward=N_backward+1
			}else {
			N_normal=N_normal+1
			}
		
		diff=interval-invl
	
}

PROCEDURE event_time() {LOCAL diff2,T
        diff2=diff
	if (number > 0) {
	   T=invl(m)
	   rnd=T
	   if(T==0 && diff2==0) { T=T+dt } 
	   
	   event = T+event + diff2    
	   
 	   N_total=N_total+1
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
		if(x == 0){  
				rand=scop_random() 
				x = rand
				net_event(t)
				event_time()  

		  		 if (event-t <= refrac+dt && event-t >= refrac) {
		      		 		event=event+dt    
		  		 }

		  		 if (on==1) {
					net_send(event - t, 1) 
		  				 }
				net_send(refrac, 2)  
		
		} else if (x!=0) {
			net_event(t)
			event_time() 
			if (on == 1) {
				net_send(event - t, 1)
					  }
		}
		
	}
	if (flag == 2) {
		
		x = 0
		
	}
}