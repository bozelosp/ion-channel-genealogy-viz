NEURON {
  SUFFIX mhw
  
  
  
  GLOBAL mode, slope_thresh
  RANGE Vthresh, Peak, tmax, Amp, dvdt 
  RANGE Halfheight, t0, t1, Width 
  RANGE Risetime, Falltime, t_rt_0, t_rt_1, t_ft_0, t_ft_1, thresh10, thresh90
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mM) = (milli/liter)

 }

PARAMETER {
  
  
  
  mode = 0 (1) 
  slope_thresh = 20 (1) 
    
}

ASSIGNED {
  v (mV)     
  Vthresh (mV) 
  VthreshSet (1) 
  Peak (mV)  
  tmax (ms)  
  Halfheight (mV) 
  t0 (ms)    
  t1 (ms)    
  Amp (mV)   
  Width (ms)    
  findwhich (1) 
  Risetime (ms)    
  Falltime (ms)    
  thresh10 (mV) 
  thresh90 (mV) 
  t_rt_0 (ms) 
  t_rt_1 (ms) 
  t_ft_0 (ms) 
  t_ft_1 (ms) 
                    
  rise10set (1) 
  rise90set (1) 
  fall90set (1)
  fall10set (1)
  
  v2 (mV)
  v3 (mV)
  t2 (mV)
  t3 (mV)
  dvdt (mV/ms)
}

INITIAL {
  if (mode==1) { 

    Vthresh = -100 (mV) 
    Peak = v
    tmax = -1 (ms) 
    Halfheight = v
    t0 = -1 (ms)
    t1 = -1 (ms)
    Width = -1 (ms)
    Risetime = -1 (ms)
    Falltime = -1 (ms)
    thresh10 = v 
    thresh90 = v
    t_rt_0 = -1 (ms)
    t_rt_1 = -1 (ms)
    t_ft_0 = -1 (ms)
    t_ft_1 = -1 (ms)
    
	v2=0 (mV)
	v3=0 (mV)
	t2=0 (ms)
	t3=0 (ms)
	dvdt=0 (mV/ms)
    VthreshSet=0 (1) 
	} else if (mode==2) { 

    Amp = Peak - Vthresh
    Halfheight = Vthresh + 0.5 * Amp 
    findwhich = 0 
    thresh10 = Vthresh + 0.1 * Amp
    thresh90 = Vthresh + 0.9 * Amp
  } else if (mode==0) {
    Vthresh = v
    Peak = v
    tmax = -1 (ms) 
    t0 = -1 (ms)
    t1 = -1 (ms)
    Width = -1 (ms)
    Risetime = -1 (ms)
    Falltime = -1 (ms)


    t_rt_0 = -1 (ms)
    t_rt_1 = -1 (ms)
    t_ft_0 = -1 (ms)
    t_ft_1 = -1 (ms)
    findwhich = 0 
  }
  rise10set = 0 (1) 
  rise90set = 0 (1) 
  fall90set = 0 (1)
  fall10set = 0 (1)
 
  }

PROCEDURE findmax() {
  if (v>Peak) {
    Peak = v
    tmax = t
  }
  if (!VthreshSet && t > 2) { 
    if (dvdt >= slope_thresh) {
      Vthresh = v
      VthreshSet = 1
     }
  }
}


PROCEDURE findx() {

  if (findwhich==0) {
    if (v > Halfheight) {
      t0 = t
      findwhich = 1
    }
  } else if (findwhich==1) {
    if (v < Halfheight) { 
      t1 = t
      Width = t1-t0
      findwhich = 2 
    }
  }

  if (rise90set) { 
    if (!fall90set) {
      if (v < thresh90) {
        t_ft_0 = t
        fall90set = 1
      }
    }
    if (!fall10set) {
      if (v < thresh10) {
        t_ft_1 = t
        fall10set = 1
        Falltime = t_ft_1 - t_ft_0
      }
    }
  }


  if (!rise10set) {
    if (v > thresh10) {
      t_rt_0 = t
      rise10set = 1
    }
  }
  if (!rise90set) {
    if (v > thresh90) {
      t_rt_1 = t
      rise90set = 1
      Risetime = t_rt_1 - t_rt_0
    }
  }
}

BREAKPOINT {
        SOLVE check METHOD after_cvode
}
PROCEDURE check() {
	v2 = v3
	v3 = v
	t2 = t3
	t3 = t
	if (t3 - t2 > 0) {
		dvdt = (v3 - v2)/(t3 - t2)
	}
}
AFTER SOLVE { 
  if (mode==1) { 
    findmax()
  } else if (mode==2) {
    findx()
  } else if (mode==0) {
    findmax() 
    findx()
  }
}