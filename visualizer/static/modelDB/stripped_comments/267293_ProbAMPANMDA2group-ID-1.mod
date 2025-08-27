NEURON {

        POINT_PROCESS ProbAMPANMDA2group
        RANGE tau_r_AMPA, tau_d_AMPA, tau_r_NMDA, tau_d_NMDA, Nsyns
        RANGE Use, u, Dep, Fac, u0, weight_factor_NMDA
        RANGE i, i_AMPA, i_NMDA, g_AMPA, g_NMDA, e, gmax
        NONSPECIFIC_CURRENT i, i_AMPA,i_NMDA
	POINTER rng
}

PARAMETER {

        tau_r_AMPA = 0.2   (ms)  
        tau_d_AMPA = 1.7    (ms)  
	tau_r_NMDA = 0.29   (ms) 
        tau_d_NMDA = 43     (ms) 
        Use = 1.0   (1)   
        Dep = 100   (ms)  
        Fac = 10   (ms)  
        e = 0     (mV)  
	mg = 1   (mM)  
        mggate
    	gmax = .001 (uS) 
    	u0 = 0 
	Nsyns = 10 
        weight_factor_NMDA = 1
}


   
VERBATIM

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#ifndef NRN_VERSION_GTEQ_8_2_0
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
extern int ifarg(int iarg);
extern int vector_capacity(void* vv);
extern void* vector_arg(int iarg);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif
ENDVERBATIM

ASSIGNED {

        v (mV)
        i (nA)
	i_AMPA (nA)
	i_NMDA (nA)
        g_AMPA (uS)
	g_NMDA (uS)
        factor_AMPA
	factor_NMDA
	rng
 	space       
}

STATE {

        A_AMPA       
        B_AMPA       
	A_NMDA       
        B_NMDA       
}

INITIAL{

        LOCAL tp_AMPA, tp_NMDA
        
	A_AMPA = 0
        B_AMPA = 0
	
	A_NMDA = 0
	B_NMDA = 0
        
	tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA) 
	tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA) 
        
	factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) 
        factor_AMPA = 1/factor_AMPA
	
	factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA) 
        factor_NMDA = 1/factor_NMDA

}

BREAKPOINT {

        SOLVE state METHOD cnexp
	mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (mg / 3.57 (mM))) 
        g_AMPA = gmax*(B_AMPA-A_AMPA) 
	g_NMDA = gmax*(B_NMDA-A_NMDA) * mggate 
        i_AMPA = g_AMPA*(v-e) 
	i_NMDA = g_NMDA*(v-e) 
	i = i_AMPA + i_NMDA
}

DERIVATIVE state{

        A_AMPA' = -A_AMPA/tau_r_AMPA
        B_AMPA' = -B_AMPA/tau_d_AMPA
	A_NMDA' = -A_NMDA/tau_r_NMDA
        B_NMDA' = -B_NMDA/tau_d_NMDA
}


NET_RECEIVE (weight, Pv, Pr, u, myInd, tsyn (ms), Pv_tmp){
	
	

        INITIAL{
                Pv=1
                u=u0
            }

        
        VERBATIM
          void** vv = (void**)(&space);
          double *x;
          int nx = vector_instance_px(*vv, &x);
          int myInd = rand()%((int)Nsyns);
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
            A_AMPA = A_AMPA + weight*factor_AMPA
            B_AMPA = B_AMPA + weight*factor_AMPA
            A_NMDA = A_NMDA + weight*weight_factor_NMDA*factor_NMDA
            B_NMDA = B_NMDA + weight*weight_factor_NMDA*factor_NMDA
            
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
                
                value = nrn_random_pick(RANDCAST _p_rng);
		        
                
                
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