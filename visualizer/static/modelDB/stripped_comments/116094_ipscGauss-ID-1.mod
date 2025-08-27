NEURON {
	POINT_PROCESS ipsc_gauss
	NONSPECIFIC_CURRENT  ip
	RANGE ip, gp
	RANGE dep, dur, delm, sigi, amp, tau, nipsc
	RANGE noise_seed
}
UNITS {
	(nA) = (nanoamp)

}

PARAMETER {
	dep  = 5 (ms)		 
	dur  = 80 (ms) <0,1e9>   
	delm = 10 (ms)		 
	sigi = 50 (ms)		 
	dt   = 0.01 (ms) 
	v (mV)			 

	amp  = 0.005    (nS)     
	tau  = 5 	(ms)     
	Erev = -70	(mV)	 


	nipsc = 5		 
	noise_seed = 1
}

PROCEDURE seed1(x) {
	set_seed(x)
}

ASSIGNED { ip (nA)
	   gp (uS)
	   indic1
	   events[20000]
	   amplitude[20000]
	   tevents[200]
	   count
}

LOCAL j,k 

INITIAL {LOCAL indic3
	ip = 0
	gp = 0
	count = 1
	indic1 = dur/dt-1
	seed1(noise_seed)
	FROM j = 0 TO indic1 {
		events[j]    = 0  
		amplitude[j] = 0  
	}

	FROM k = 0 TO nipsc-1 {
		indic3 = -1
		WHILE (indic3<0 || indic3>dur/(dt)) {
		      tevents[k] = normrand(delm,sigi)
		      indic3 = floor((tevents[k]-dep)/(dt)) 
		}
		amplitude[indic3]   = scop_random() 
		
	        events[indic3]      = events[indic3] + amplitude[indic3]  
	}
}

BREAKPOINT { LOCAL indic2
	if (t<dep+dur && t >= dep) {
		indic2 = floor((t-dep)/(dt)) 
		gp = (-gp/tau)*(dt/2) + amp * (amplitude[indic2]/2) * events[indic2] + gp 
	
		ip = gp * (v - Erev)
		
		
		if (ip<0 || ip>100) {
			ip = 0			
		}
		if (ip>10) {
			ip = 10
		}
	}
	else {
		gp = 0
	}
}