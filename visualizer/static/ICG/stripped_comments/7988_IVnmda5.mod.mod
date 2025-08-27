INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS IVNMDA5
	POINTER C
	RANGE C0, C1, C2, D, O, B
	RANGE g, gmax, rb, gluc, w, timer, llabel
	GLOBAL Erev, mg, Rb, Ru, Rd, Rr, Ro, Rc
	GLOBAL vmin, vmax

      NONSPECIFIC_CURRENT inmda
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= 0    (mV)	
	gmax	= 500  (pS)	
	mg	= 0    (mM)	
	vmin = -120	(mV)
	vmax = 100	(mV)

	gluc = 0.00003 (mM)
	timer = 110 (ms)
	llabel = 0
	
	
	Rb	= 5e-3    (/uM /ms)	
	Ru	= 12.9e-3  (/ms)	
	Rd	= 8.4e-3   (/ms)	
	Rr	= 6.8e-3   (/ms)	
	Ro	= 46.5e-3   (/ms)	
	Rc	= 73.8e-3   (/ms)	
}




ASSIGNED {
	v		(mV)		
        inmda               (nA)            
	g 		(pS)		
	C 		(mM)		

	rb		(/ms)    
	w
}

STATE {
	
	C0		
	C1		
	C2		
	D		
	O		

	B		
}

INITIAL {
	rates(v)
	C0 = 1
	llabel = 0.0
}

BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse

	g = gmax * O * B
      inmda = (1e-6) * g * (v - Erev)



	
}

KINETIC kstates {
	
	rb = Rb * (1e3) * C 

        ~ C0 <-> C1     (2*rb,Ru)
        ~ C1 <-> C2     (rb,2*Ru)
	~ C2 <-> D	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)

	CONSERVE C0+C1+C2+D+O = 1
}

PROCEDURE rates(v(mV)) {
	TABLE B
	DEPEND mg
	FROM vmin TO vmax WITH 200

	

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

VERBATIM
double peak = 0.0;
double writeIV(double(in),double(cc), double(timerr), double(label)) {

			 
                      FILE *fp;


			    if (peak*1e10 >= in*1e10) {
					peak = in;

					if(label == 1.0) {
						label = 0.0;
					}
					if(label == 2.0) {
						label = 0.0;
					}
					if(label == 3.0) {
						peak = 0;
					}
								
			    }

			    if ((peak*1e10 < in*1e10) && (label != 3.0)) {
  			         
				   label = label + 1.0;
				   if(label == 1.0){
					peak = peak;
				   } else if (label == 3.0) {
    					   fp = fopen("IV.txt", "a");

	                      	   if(fp == NULL)  {
#if 0
						      printf("Can't open IV.txt file\n");
						      exit(1);
#else
							hoc_execerror("Can't open IV.txt file", (char*)0);
#endif
	                           }
					
					   printf("peak value is %e\n", peak);	
					   fprintf(fp, "%e	%e \n", cc, peak);
					   fclose(fp);
					   peak = 0.0;
				   }				
			   }
			   if ((t+dt/2) > timerr ) {
				printf("baby one");
		            if (label == 0.0) {
  			         fp = fopen("IV.txt", "a");
                      	   if(fp == NULL)  {
#if 0
					      printf("Can't open IV.txt file\n");
					      exit(1);
#else
							hoc_execerror("Can't open IV.txt file", (char*)0);
#endif
	                     }
				   fprintf(fp, "%e %e \n", cc, in);
					   
                           fclose(fp);
				   label = 3.0;
				   peak = 0.0;		
				}	
					  
			   }	
			                			
				return(label);		
}
ENDVERBATIM