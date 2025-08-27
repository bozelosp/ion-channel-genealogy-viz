NEURON {

        POINT_PROCESS ProbUDFsyn2groupdet  
        RANGE tau_r, tau_d, Nsyns, Nevents, eventCounter
        RANGE Use, u, Dep, Fac, u0
        RANGE i, g, e, gmax
        NONSPECIFIC_CURRENT i
	POINTER rng
}

PARAMETER {

        tau_r = 0.2   (ms)  
        tau_d = 1.7    (ms)  
        Use = 1.0   (1)   
        Dep = 100   (ms)  
        Fac = 10   (ms)  
        e = 0     (mV)  
    	gmax = .001 (uS) 
    	u0 = 0 
        Nsyns = 10 
        Nevents = 0 
}


   
VERBATIM

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);

extern int ifarg(int iarg);
extern int vector_capacity(void* vv);
extern void* vector_arg(int iarg);

ENDVERBATIM
  

ASSIGNED {

        v (mV)
        i (nA)
	g (uS)
        factor
	rng
	weight_NMDA
        space        
        space2       
        eventCounter 
}

STATE {
        A       
        B       
}

INITIAL{

  LOCAL tp
        
	A = 0
  B = 0
	
        
	tp = (tau_r*tau_d)/(tau_d-tau_r)*log(tau_d/tau_r) 
	      
	factor = -exp(-tp/tau_r)+exp(-tp/tau_d) 
        factor = 1/factor 
        eventCounter = 0

}

BREAKPOINT {

        SOLVE state METHOD cnexp
        g = gmax*(B-A) 
        i = g*(v-e) 
}

DERIVATIVE state{

        A' = -A/tau_r
        B' = -B/tau_d
}


NET_RECEIVE (weight, Pv, Pr, u, myInd, tsyn (ms), Pv_tmp){
	
        INITIAL{
                Pv=1
                u=u0
                tsyn=t
		eventCounter=0
            }

        
        VERBATIM
          void** vv = (void**)(&space);
          void** vv2 = (void**)(&space2);
          double *x;
          int nx = vector_instance_px(*vv, &x);
          double *x2;
          int nx2 = vector_instance_px(*vv2, &x2);
          int myInd = 0;
          if (eventCounter < nx2)
            myInd = x2[(int)eventCounter];
          else printf("eventCounter >= nx2! t = %lf\n",t);
          eventCounter++;
          _args[4] = myInd;
          _args[5] = x[myInd];                
          _args[1] = x[myInd+(int)Nsyns];     
          _args[3] = x[myInd+2*((int)Nsyns)]; 
        ENDVERBATIM
        
        

        
        if (Fac > 0) {
          u = u*exp(-(t - tsyn)/Fac) 
        } else {
          u = Use  
        } 
        if(Fac > 0){
          u = u + Use*(1-u) 
        }    

        
        Pv_tmp  = 1 - (1-Pv) * exp(-(t-tsyn)/Dep) 
                                                  
        Pr  = u * Pv_tmp                          
        Pv_tmp  = Pv_tmp - u * Pv_tmp             
        
        
                
	if (erand() < Pr){
          tsyn = t
	  Pv = Pv_tmp
          A = A + weight*factor
          B = B + weight*factor
          
          
          
        } else {
          
          
          
        }
        VERBATIM
          x[myInd] = _args[5];
          x[myInd+(int)Nsyns] = _args[1];
          x[myInd+2*((int)Nsyns)] = _args[3];
        ENDVERBATIM
}

PROCEDURE setRNG() {
VERBATIM
    {
        
        void** pv = (void**)(&_p_rng);
        if( ifarg(1)) {
            *pv = nrn_random_arg(1);
        } else {
            *pv = (void*)0;
        }
    }
ENDVERBATIM
}

FUNCTION erand() {
VERBATIM
	    
        double value;
        if (_p_rng) {
                
                value = nrn_random_pick(_p_rng);
		        
                
                
                return value;
        }else{
ENDVERBATIM
                
                
                
                
                erand = exprand(1)
VERBATIM
        }
ENDVERBATIM
        
                       
                       
}

PROCEDURE setVec() {    
                        
  VERBATIM
  void** vv;
  vv = (void**)(&space);
  *vv = (void*)0;
  if (ifarg(1)) {
    *vv = vector_arg(1);
    Nsyns = vector_capacity(*vv)/3;
  }
  ENDVERBATIM
}

PROCEDURE printVec() { 
VERBATIM
    void** vv = (void**)(&space);
    double *x;
    int nx = vector_instance_px(*vv, &x);
    int i1;
    for (i1=0; i1<Nsyns;i1++) {
      printf("tsyns[%i] = %g, Pv[%i] = %g, u[%i] = %g\n", i1, x[i1], i1, x[i1+(nx/3)], i1, x[i1+2*(nx/3)]);
    }
ENDVERBATIM
}

PROCEDURE setVec2() {    
  VERBATIM
  void** vv;
  vv = (void**)(&space2);
  *vv = (void*)0;
  if (ifarg(1)) {
    *vv = vector_arg(1);
    Nevents = vector_capacity(*vv);
  }
  ENDVERBATIM
}

PROCEDURE printVec2() { 
VERBATIM
    void** vv = (void**)(&space2);
    double *x;
    int nx = vector_instance_px(*vv, &x);
    int i1;
    for (i1=0; i1<Nevents;i1++) {
      printf("%g ", x[i1]);
      if (i1%100==99)
        printf("\n");
    }
ENDVERBATIM
}