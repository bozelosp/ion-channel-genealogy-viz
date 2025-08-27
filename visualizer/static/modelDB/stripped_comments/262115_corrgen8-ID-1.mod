DEFINE MAXCHANNELS 25000		

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON	{ 
	POINT_PROCESS corrGen8
	RANGE N, freq, correl, refract, min_val, max_val
	RANGE on, latency, shutoff
	RANGE x, ns, ls, sync
	RANGE spont, prob, N2
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt	(ms)

	N	= 100		
	N2	= 1		
	freq	= 40 (/s)	
	correl	= 0		
	refract = 1 (ms)	
	min_val	= -70 (mV)	
	max_val	= 40 (mV)	
	on	= 1		
	latency	= 0 (ms)	
	shutoff	= 1e6 (ms)	

}

ASSIGNED {
	x[MAXCHANNELS]	(mV)	
	R[MAXCHANNELS]		
	ns[MAXCHANNELS]		
	ls[MAXCHANNELS]		
	sync[MAXCHANNELS]	
	spont			
	prob			
}
	
INITIAL { LOCAL i
	spont = (0.001) * freq * dt	
	if(spont > 0.5) {
	   VERBATIM
	   printf("\n\nERROR in correlated random generator\n");
	   printf("firing probability is too high: spont=%g\n",(float)spont);
	   printf("decrease integration step or mean firing frequency\n");
	   exit(-1);
	   ENDVERBATIM
	}
	FROM i=0 TO N-1 {
	  ns[i] = 0
	  ls[i] = -10000
	}
	N2 = N + sqrt(correl) * (1-N)
	go()
}

BREAKPOINT {
	SOLVE go
}


UNITSOFF

PROCEDURE go() { LOCAL i, j, sum



   FROM i=0 TO N-1 {
	x[i] = min_val
   }

   if( (on==1) && (t>=latency) && (t<=shutoff) )  {



     if(correl==0) {		

	FROM i=0 TO N-1 {
	   if(get_random(1) <= spont) { 	
		x[i] = max_val			
	   } else {
		x[i] = min_val
	   }
	}

     } else {			

	FROM i=0 TO N2-1 {		
	   if(get_random(1) <= spont) { 	
		R[i] = max_val			
	   } else {
		R[i] = min_val
	   }
	}

	FROM i=0 TO N-1 {		
	   j = get_random(N2)			
	   x[i] = R[j]				
	}

     }



     FROM i=0 TO N-1 {
	   if((t-ls[i]) <= refract) {		
		x[i] = min_val
	   }
	   if(x[i] == max_val) {		
		   ns[i] = ns[i] + 1			
		   ls[i] = t				
	   }
     }


   } else if(on==2) {
	FROM i=0 TO N-1 {		
	   x[i] = max_val				
	   ns[i] = ns[i] + 1				
	   ls[i] = t					
	}
	on = 1				

   } else if(on==3) {
	FROM i=0 TO N-1 {		
	   if(sync[i]) {
		x[i] = max_val				
		ns[i] = ns[i] + 1			
		ls[i] = t				
	   }
	}
	on = 1				
   }
}
UNITSON


FUNCTION get_random(maxval) {			
	get_random = maxval * random() / (2^31)
}


PROCEDURE new_seed(seed) {		
	srandom(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}

PROCEDURE printvec() { LOCAL i		
   VERBATIM 
   {
	int i;
	printf("i\tx\tns\tls\n");
	for(i=0; i<N; i++) {
	  printf("%d\t%g\t%g\t%g\n",i,(float)x[i],(float)ns[i],(float)ls[i]);
	}
   } 
   ENDVERBATIM
}