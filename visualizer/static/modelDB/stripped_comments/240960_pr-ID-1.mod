NEURON {
	POINT_PROCESS Pr
	RANGE P, P0, random, f, tau_F, d1, tau_D1, F, D1, tlast
    THREADSAFE
    POINTER randObjPtr
}

PARAMETER {
    
    P0 = 0.200          (1)     < 0, 1 >        
    f = 1.769           (1)     < 0, 1e9 >      
    tau_F = 67.351      (ms)    < 1e-9, 1e9 >   
    d1 = 0.878          (1)     < 0, 1 >        
    tau_D1 = 92.918     (ms)    < 1e-9, 1e9 >   
    
    
}

ASSIGNED {
	P				        
    randObjPtr              
    random                  
    F                       
    D1                      
    
    tlast (ms)              
}

INITIAL {
	P = P0
    random = 1
    F = 1
    D1 = 1
    
}

NET_RECEIVE(weight) {
    INITIAL {
        tlast = t
    }
    F = 1 + (F-1)*exp(-(t - tlast)/tau_F)
    D1 = 1 - (1-D1)*exp(-(t - tlast)/tau_D1)


    if (P0*F*D1 > 1) {
        P = 1
    } else {

        P = P0*F*D1
    }
    random = randGen()
    if (random <= P) {
        net_event(t)
    }
    tlast = t
    F = F + f
    D1 = D1 * d1

}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION randGen() {
VERBATIM
   if (_p_randObjPtr) {
      
      _lrandGen = nrn_random_pick(_p_randObjPtr);
   }else{
      hoc_execerror("Random object ref not set correctly for randObjPtr"," only via hoc Random");
   }
ENDVERBATIM
}

PROCEDURE setRandObjRef() {
VERBATIM
   void** pv4 = (void**)(&_p_randObjPtr);
   if (ifarg(1)) {
      *pv4 = nrn_random_arg(1);
   }else{
      *pv4 = (void*)0;
   }
ENDVERBATIM
}