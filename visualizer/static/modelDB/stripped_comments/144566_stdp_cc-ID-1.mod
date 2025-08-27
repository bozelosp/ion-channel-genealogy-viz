NEURON {
       POINT_PROCESS STDPSynCC
       RANGE tau, e, i 
       RANGE A_m, A_p, tau_y, tetam,tetap,tau_r,tau_0 
       RANGE t_last_pre 
       RANGE u_m1, u_m2,r,g_update,gbar
       RANGE delay_steps,delay_array_u_m1, delay_array_u_m2, pointless_counter 
       RANGE um2s,um1s
       NONSPECIFIC_CURRENT i
}

DEFINE MAX_ARRAY 2000 

UNITS {
      (nA) = (nanoamp)
      (mV) = (millivolt)
      (uS) = (microsiemens)
}

PARAMETER {
	  
	  tau = 0.1 (ms) <1e-9,1e9>
	  e = 0	(mV)
	  gbar = 0.05

	  
	  tau_0 = 10 
	  tau_r =15 
	  tau_y = 114 
	  A_m = 0.00001 
	  A_p = 0.00012
	  tetam = -64.9 
	  tetap = -35 
	  delay_steps = 50 

	  
	  um1s
	  um2s

	  
	  t_last_pre = -1
	  g_update
	  pointless_counter = 0
}

ASSIGNED {
	 v (mV)
	 i (nA)
	 delay_array_u_m1[MAX_ARRAY]
	 delay_array_u_m2[MAX_ARRAY]
	 delay_array_v[MAX_ARRAY]
}

STATE {
      g (uS)
      u_m1
      u_m2
      r
}

INITIAL {
	g=0
	u_m1 = v 
	u_m2 = v 
	r = 0
	FROM i=0 TO MAX_ARRAY {
	     delay_array_v[i] = v
	     delay_array_u_m1[i] = 0
	     delay_array_u_m2[i] = 0
	}
}

BREAKPOINT {
	   SOLVE state METHOD cnexp
	   i = g*(v - e)
}

DERIVATIVE state { 
	   LOCAL x,u_sig,u_m1_sig,u_m2_sig,set_loc,retrieve_loc 

	  
	  if( (v - tetap) >  0) {
	      u_sig = v - tetap 
	      }
	  else { u_sig = 0 }

	  if( (delay_array_u_m1[retrieve_loc] - tetam) > 0) {
	      u_m1_sig = delay_array_u_m1[retrieve_loc] - tetam 
	      um1s = u_m1_sig
	  }
	  else { 
	       u_m1_sig = 0 
	       um1s = 0
	  }
	  if( (delay_array_u_m2[retrieve_loc] - tetam) > 0 ) {
	      u_m2_sig = delay_array_u_m2[retrieve_loc] - tetam 
	      um2s = u_m2_sig
	  }
	  else { 
	       u_m2_sig = 0 
	       um2s = 0
	  }

	  if(t_last_pre == 1) {
	  	x = 1
	    	t_last_pre = 0.5
	    	
	  } else {
	    	x = 0
	  }

	  g' = -g/tau
	  u_m1' = (v-u_m1)/tau_0 
	  u_m2' = (v-u_m2)/tau_y
	  r' = (x-r)/ tau_r
	  g_update = - A_m*x*u_m1_sig + A_p*u_sig*r*u_m2_sig 
	  gbar = gbar + g_update

	  
	  set_loc = fmod(pointless_counter,delay_steps)
	  retrieve_loc = fmod(pointless_counter+delay_steps+1,delay_steps)
	  delay_array_u_m1[set_loc] = u_m1
	  delay_array_u_m2[set_loc] = u_m2
	  pointless_counter = pointless_counter + 1 
}


NET_RECEIVE(weight (uS)) {
		   g = gbar 
		   t_last_pre = 1
		   printf("received weight=%f, at t=%f\n", weight,t )
}