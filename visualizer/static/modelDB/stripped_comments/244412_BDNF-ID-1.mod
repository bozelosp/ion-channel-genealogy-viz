NEURON {
    POINT_PROCESS BDNF
    USEION ca READ cai	
    
    
    RANGE max_BDNF_rel_delay, theta_cai_BDNF, max_cai_BDNF, BDNF_prel, fused_vesicles, duration_BDNF_release
    RANGE proBDNF_uptake,PC_uptake,mBDNF_uptake, TrkB, proBDNF_fraction
    RANGE tau_LTP14, theta_gAMPA, sigma_gAMPA, alpha_gAMPA, shift_gAMPA, v_BDNF, is, intracell_signaling
    POINTER randObjPtr
    POINTER gAMPA
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {
    cai (mM)
    
    dt (ms)
       
    
    
    theta_cai_BDNF = 0.11 (mM)
    max_cai_BDNF = 0.13 (mM)
    
    max_BDNF_rel_delay = 400e3 (ms)
    duration_BDNF_release = 1800e3 
    proBDNF_uptake = 0.00001 (mM/ms)    
    mBDNF_uptake = 0.00001 (mM/ms)    
    PC_uptake = 0.00001 (mM/ms)    
    
    
    cai_integration_time_step = 1 (ms)
    
    alpha_dt_cai = 1.2 (1/mM 1/ms)
    
    
    n_vesicles = 300 (1)
    v_BDNF = 0.002 (mM)
    proBDNF_fraction = 0.7
    v_PC = 0.002 (mM)
    
    
    tau_cleave = 10000.0 (ms/mM) 
    
    
    
    tau_LTP14 = 180e3 
    theta_TrkB = 0.0002 (mM)
    sigma_TrkB = 0.00001 (mM)
    
    gbar_AMPA = 1 (nS)
    scale_AMPA = 100 (1)
    alpha_gAMPA = 0.5
    theta_gAMPA = 0.5
    sigma_gAMPA = 0.1
    shift_gAMPA = 0    
    
    
    
    
    is_decay = 0.1e-3 (/ms) 
    
    
}

ASSIGNED {
    BDNF_prel (1)
    randObjPtr
    cai_th_crossed (1)
    gAMPA (nS)
    gAMPA_g
    alpha_LTP14 (mM/ms)
    mBDNF_fraction
    fusion_delay (ms)
    v_prel_norm (1)
    n_avail_vesicles (1)
    cai_factor (1)
    
}

STATE {
    
    
    fused_vesicles (1)
    proBDNF (mM)
    mBDNF (mM)
    ppBDNF (mM) 
    PC (mM)  
    proBDNF_PC (mM)
    TrkB (mM)
    proBDNF_removed (mM)
    mBDNF_removed (mM)
    PC_removed (mM)
    is (1)
    intracell_decay (1)
    intracell_signaling (mM)
}

INITIAL {    
    fused_vesicles = 0
    cai_th_crossed = 0
        
    proBDNF = 0 
    mBDNF = 0 
    ppBDNF = 0
    PC = 0
    proBDNF_removed = 0 
    mBDNF_removed = 0 
    ppBDNF = 0
    PC_removed = 0
    is = 0    
    intracell_decay = 0
    intracell_signaling = 0
    mBDNF_fraction = 1 - proBDNF_fraction

    alpha_LTP14 = 1/tau_LTP14
    
    TrkB = 0 (mM)
    
    n_avail_vesicles = n_vesicles
}

BREAKPOINT {
    SOLVE kstates METHOD sparse
    
    
}

KINETIC kstates {
    
    
    
    
    
    
    
    
    
    
    ~ is <-> intracell_decay (is_decay,0)
    
    
    
    
    
    ~ proBDNF << (v_BDNF * proBDNF_fraction * fused_vesicles / duration_BDNF_release)
    ~ mBDNF << (v_BDNF * mBDNF_fraction * fused_vesicles / duration_BDNF_release)
    ~ PC << (v_PC * fused_vesicles / duration_BDNF_release)
    
    
    ~ proBDNF + PC <-> proBDNF_PC (1/tau_cleave,0)  
    ~ proBDNF_PC <-> ppBDNF + mBDNF + PC (1/tau_cleave,0)  
    
    ~ proBDNF <-> proBDNF_removed (proBDNF_uptake, 0)
    ~ mBDNF <-> mBDNF_removed (mBDNF_uptake, 0)
    ~ PC <-> PC_removed (PC_uptake, 0)
    ~ fused_vesicles <-> fused_vesicles (0,0)

    
    TrkB = mBDNF * sigh(mBDNF, theta_TrkB, sigma_TrkB)
    ~ TrkB <-> intracell_signaling (alpha_LTP14, 0)
    
    
    
    gAMPA = 1 + alpha_gAMPA * sigh(intracell_signaling, theta_gAMPA, sigma_gAMPA) + shift_gAMPA
    
}


FUNCTION sigh(x (mM), theta (mM), sigma (mM)) {
    
    
    
    
    
    sigh = 1 / (1 + exp((theta - x) / sigma))
}

VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif
ENDVERBATIM

FUNCTION randGen() {
VERBATIM
   if (_p_randObjPtr) {
      
      _lrandGen = nrn_random_pick(RANDCAST _p_randObjPtr);
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

FUNCTION min(x,y) { if (x<=y){ min = x }else{ min = y } }

NET_RECEIVE (weight (1)) { LOCAL prel, is_effect
    if ((flag == 0) && (cai_th_crossed == 0) ) { 
	cai_th_crossed = 1
	
	is = is + 0.1
	
	
	
	
	if (n_avail_vesicles > 0) {
	    net_send(cai_integration_time_step,2) 
	}
    }
    
    if (flag == 2 && n_avail_vesicles > 0) {
    	
	cai_factor = (cai - theta_cai_BDNF) / (max_cai_BDNF - theta_cai_BDNF)
	if (is > 0.15) {
	    is_effect = 1
	} else {
	    is_effect = 0
	}
	prel = alpha_dt_cai * cai_integration_time_step * min(1,cai_factor) * is_effect
	
	
    	if ( randGen() < prel ) {
	    BDNF_prel = prel
    	    
    	    
	    
	    fusion_delay = max_BDNF_rel_delay * (1 - min(1,cai_factor)) * randGen()
    	    net_send(fusion_delay, 3)
	    n_avail_vesicles = n_avail_vesicles - 1
    	    
    	    
    	    
    	}
	
    	if (cai > theta_cai_BDNF) {
	    if (n_avail_vesicles > 0) {
    		net_send(cai_integration_time_step,2) 
	    }
	} else {
	    cai_th_crossed = 0
	}
    }
    
    if (flag == 3) { 
	
	fused_vesicles = fused_vesicles + 1
	
	net_send(duration_BDNF_release,4)
    }
    
    if (flag == 4) { 
	
	fused_vesicles = fused_vesicles - 1
    }
}