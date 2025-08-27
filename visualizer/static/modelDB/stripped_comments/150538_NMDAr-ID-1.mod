INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
     POINT_PROCESS NMDA_TESTED
     RANGE state_C0, state_C1, state_C2, state_D, state_O, state_B,  
state_DB, state_C2B, state_C1B, state_CB
     RANGE g, nchan, gamma, rb, mg_on, mg_off, C, conc, del, dur
     GLOBAL Erev, mg, kon, koff, kd, kr, beta, alpha
     GLOBAL beta_mg, alpha_mg, koff_mg, rb_mg, kon_mg, kr_mg, kd_mg
     GLOBAL vmin, vmax
     NONSPECIFIC_CURRENT i
}

UNITS {
     (nA) = (nanoamp)
     (mV) = (millivolt)
     (uS) = (microsiemens)
     (umho) = (micromho)
     (mM) = (milli/liter)
     (uM) = (micro/liter)
}

PARAMETER {

     Erev    = -3    (mV) 
     nchan  = 10          
     gamma = 0.00005 (uS) 
     mg  = 1     (mM)     
     vmin = -180 (mV)
     vmax = 100  (mV)
     normfactor = .06948171



     
     kon  = 5e-3    (/uM /ms)    
     koff  = 82e-3  (/ms)      
     beta  = 46.5e-3   (/ms)     
     alpha  = 91.6e-3   (/ms)    
     kr  = 1.8e-3   (/ms)        
     kd  = 8.4e-3   (/ms)        

     kon_mg  = 5e-3    (/uM /ms) 
     koff_mg  = 82e-3  (/ms)   
     beta_mg  = 41.85e-3   (/ms)  
     alpha_mg  = 229e-3   (/ms) 
     kr_mg  = 1.8e-3   (/ms)     
     kd_mg  = 8.4e-3   (/ms)     
     del = 100 (ms)
     dur = 10 (ms)
     conc = 0       (uM)   
}

ASSIGNED {
     v       (mV)        
     i       (nA)        
     g       (uS)        
     rb      (/ms)       
     rb_mg   (/ms)   
     mg_on   (/mM /ms)       
     mg_off  (/ms)       
     C       (uM)       
}

STATE {
     
     state_C0       
     state_C1       
     state_C2       
     state_D        
     state_O        
     state_B        
     state_DB      
     state_C2B     
     state_C1B     
     state_CB      
}

INITIAL {
SOLVE kstates STEADYSTATE sparse
}

BREAKPOINT {
     transmitter()
     SOLVE kstates METHOD sparse

     g = nchan * state_O * gamma 
     i = g * (v - Erev) 
}

KINETIC kstates {
	rates(v)
     rb = kon * C
     rb_mg = kon_mg * C

     ~ state_C0  <-> state_C1   (2*rb,koff)
     ~ state_C1  <-> state_C2   (rb,2*koff)
     ~ state_C2  <-> state_D    (kd,kr)
     ~ state_C2  <-> state_O    (beta,alpha)
     ~ state_O   <-> state_B    (mg_on,mg_off)
     ~ state_C2B <-> state_DB   (kd_mg,kr_mg)
     ~ state_B   <-> state_C2B  (alpha_mg,beta_mg)
     ~ state_C2B <-> state_C1B  (2*koff_mg,rb_mg)
     ~ state_C1B <-> state_CB   (koff_mg,2*rb_mg)

     CONSERVE  
state_C0+state_C1+state_C2+state_D+state_O+state_B+state_DB+state_C2B+state_C1B+state_CB = 1
}

PROCEDURE rates(v(mV)) {

     

     mg_on = 610*exp(-v/17)*(1/1000)
     mg_off = 5.4*exp(v/47)
}












NET_RECEIVE(w (uM)) {

  
           
       del=t
       conc=w

}

PROCEDURE transmitter() {
if ((conc!=0)&&(t>=del)&&(t<=del+dur)) {
   C=conc
} else {
   C=0
   conc=0
}
}