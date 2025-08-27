UNITS {
    (molar) = (/liter)
    (mM) = (millimolar)
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (uS) = (micromho)
}

NEURON {
    SUFFIX nmda
    POINT_PROCESS nmda
    NONSPECIFIC_CURRENT i
    RANGE onset, gmax, e, i, g, genv, q,
          Mg, taup, pinf,
          alf, alfA, alfslope,
          bet, betA, betslope,
          tau_on, tau_on0, tau_onslope,
          tau_off1, tau_off1_0, tau_off1slope,
          tau_off2, tau_off2_0, tau_off2slope,
          f_fast, f0, fslope,
          q10, Mg_time_factor
}

PARAMETER {
    celsius= 32 
    q10 = 3 ()
    Mg_time_factor = 1 ()
    dt (ms)
    onset=0 (ms)
    e=0	(mV)
    v	(mV)
    Mg= 1.8
    gmax = 0 (umho)  
    alfslope = 47 (mV)
    alfA = 5.4 (/ms)
    betslope = 17 (mV)
    betA = 0.61(/mM-ms)
    tau_on0 = 2.915 (ms)
    tau_onslope = -0.004125 (ms/mV)
    tau_off1_0 = 61.5 (ms)
    tau_off1slope = 0.5625 (ms/mV)
    tau_off2_0 = 352.5 (ms)
    tau_off2slope = 5.7375 (ms/mV)
    f0 = 0.515  
    fslope = -0.003125 (/mV)
}

ASSIGNED {
i (nA)  g (uS) genv (uS) pinf alf bet
taup (ms) tau_on (ms) tau_off1 (ms) tau_off2 (ms) f_fast
}

STATE {
 p q (nanocoulombs)
}

UNITSOFF

INITIAL {
   nmda_rates(v)
   nmda_taus(v,1.1*onset)
   p = pinf
   q = 0  
}

BREAKPOINT {
    LOCAL trel
    SOLVE nmda_states METHOD cnexp
    if (t>onset) {
        trel= t - onset
    	genv = gmax*( -            exp(-trel/tau_on)
                   +     f_fast*exp(-trel/tau_off1)
                   + (1-f_fast)*exp(-trel/tau_off2)
                 ) 
        g=genv*p  
    } else {
        genv=0
        g=0
    }
    i = g*(v - e)
}

DERIVATIVE nmda_states {
      nmda_rates(v) 
      nmda_taus(v,t)
      p' = (pinf - p)/taup
      q' = i*(1e-3)  
}

PROCEDURE nmda_taus(v,t) { 
        LOCAL temp_factor
        temp_factor = q10^((celsius - 28.50)/10)
        if (t>onset) {
           f_fast= f0 + fslope*v
           if (f_fast>1) {
             f_fast=1
           }
           if (f_fast<0) {
             f_fast = 0
           }
           tau_on  =(tau_on0 + tau_onslope*v)/temp_factor
           tau_off1=(tau_off1_0 + tau_off1slope*v)/temp_factor
           tau_off2=(tau_off2_0 + tau_off2slope*v)/temp_factor
           if (tau_off1<tau_off1_0/10) {
              tau_off1=tau_off1_0/10
           }
           if (tau_off2<tau_off1) {
              tau_off2=tau_off1
           }
        }
}

PROCEDURE nmda_rates(v) { 
        LOCAL  temperature_factor 
        TABLE pinf, taup, alf, bet
          DEPEND q10, celsius,
                 alfA, betA, Mg, alfslope, betslope
          FROM -100 TO 100 WITH 200
        temperature_factor = q10^((celsius - 20)/10)
        alf = temperature_factor*alfA*exp(v/alfslope)
        bet = temperature_factor*betA*Mg*exp(-v/betslope)
        
	taup = Mg_time_factor/(alf + bet)
        
        pinf = alf/(alf + bet)
}

UNITSON