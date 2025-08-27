NEURON {
                                    SUFFIX kca
                                    USEION k READ ek WRITE ik
                                    USEION ca READ cai
                                    RANGE cac, gk, gbar, m_inf, tau_m,ik, cainit
                                    GLOBAL beta
                                    
                            }


                            UNITS {
                                    (mA) = (milliamp)
                                    (mV) = (millivolt)
                                    (molar) = (1/liter)
                                    (mM) = (millimolar)
                            }


                            PARAMETER {
									cai (mM) 
                                    v               (mV)
                                    celsius = 36    (degC)
                                    ek      = -80   (mV)
                                    cainit     = 100e-6 (mM)          
                                    gbar    = 0.01   (mho/cm2)
                                    beta    = 0.03   (1/ms)          
                                    cac     = 0.025  (mM)            
       				    taumin  = 5    (ms)            
                                    gk
                                  }


                            STATE {m}        

                            ASSIGNED {       
                                    ik      (mA/cm2)
                                    m_inf
                                    tau_m   (ms)
                                    tadj
                                    
                            }
                            BREAKPOINT { 
                                    SOLVE states METHOD derivimplicit
                                    gk = gbar*m*m*m     
                                    ik = gk*(v - ek)    
                            }

                            DERIVATIVE states { 
                                    evaluate_fct(v,cai)
                                   m' = (m_inf - m) / tau_m
                            }

                            UNITSOFF
                            INITIAL {
                            
                            
                            
                            
                                    tadj = 3 ^ ((celsius-22.0)/10) 
                                    evaluate_fct(v,cainit)
                                    m = m_inf
                            }

                            PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car
                                    car = (cai/cac)^2
                                    m_inf = car / ( 1 + car )      
                                    tau_m =  1 / beta / (1 + car) / tadj
                                    if(tau_m < taumin) { tau_m = taumin }   
                            }
                            UNITSON