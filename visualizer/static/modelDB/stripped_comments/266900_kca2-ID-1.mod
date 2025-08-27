NEURON {
                                    SUFFIX kca
									POINTER stim_i
                                    USEION k READ ek WRITE ik
                                    USEION ca READ cai
                                    RANGE flag, curr, gk, gbar, m_inf, tau_m,ik,ek2,vrun,count,vvrun,taun2,vrun2,delta2, stim_moltK
                                    GLOBAL  beta, cac,alpha
                            }


                            UNITS {
                                    (mA) = (milliamp)
                                    (mV) = (millivolt)
                                    (molar) = (1/liter)
                                    (mM) = (millimolar)
                            }


                            PARAMETER { 
							      curr
								  
                                    v               (mV)
                                    celsius = 36    (degC)
                                    ek      = -80   (mV)
									 ek2      = -80   (mV)
                                    cai     = 2.4e-5 (mM)           
                                    gbar    = 0.01   (mho/cm2)
                                    beta    = 0.03 
                                    cac     = 0.035  (mM)            
       				               taumin  = 0.5    (ms)            
                                    gk
									count=1
									vrun (mV)
									delta=0
									vinit=-76.2
									alpha=1.06
									 timestep=1000
									vrun2
									v0
									dv0
									ddv
									flag=0
									FK = 2
									PK = 1
									BK = 2.11
									CK = 48
									stim_moltK=1
																		
                                  }


                            STATE {m}        

                            ASSIGNED {       
                                    ik      (mA/cm2)
                                    m_inf
                                    tau_m   (ms)
                                    tadj
									vvrun
									taun2
									stim_i
                            }
                            BREAKPOINT { 
                                    SOLVE states METHOD cnexp
                                    gk = gbar*m*m*m     
									ek2=ek+vvrun*alpha
									ik = gk*(v - ek2)    
                            }

                            DERIVATIVE states { 
                                    evaluate_fct(v,cai)
                                    m' = (m_inf - m) / tau_m
                            }

                            UNITSOFF
                            INITIAL {
                            
                            
                            
                            
							        vrun=0
									vvrun=vrun
                                    tadj = 3 ^ ((celsius-22.0)/10) 
                                    evaluate_fct(v,cai)
                                    m = m_inf
                            }
		
		BEFORE STEP { LOCAL i
       	
		  if(stim_i==0 && flag==0){ 
		  vrun=0
		  vvrun=0
		  
	    }else{
		 flag=1
		             		  
		delta=v-vinit
		if (count<timestep+1){
		   vrun= (delta-vrun)*(FK/(count+1))+vrun
	       vrun2=vrun 
		 }else{

		vrun2= (delta)*(FK/(timestep+1))+vrun2*pow((1-FK/(timestep+1)),PK)
			
			}
		
	   vvrun=(BK*vrun2/(1+vrun2/CK))
	    
		count=count+1   
        }						
		 
	
}
   
   
                            PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car,i
                                    car = (cai/cac)^4
                                    m_inf = car / ( 1 + car )      
                                    tau_m =  1 / beta / (1 + car) / tadj
                                    if(tau_m < taumin) { tau_m = taumin }   
									   
									    
		
				
									
									
								}
                            UNITSON