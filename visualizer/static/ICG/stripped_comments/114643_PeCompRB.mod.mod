VERBATIM
#define kernel_filename	"kernel.txt"
ENDVERBATIM



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS eCompRB
	POINTER vcell,i_out
	RANGE active,active_pulses
	RANGE vr, vr10, vr50, ue, ue10, ue50
	RANGE kernel
	RANGE ksize
	RANGE imax
	RANGE t_start

RANGE g_e,g_i,g_eout,g_iout,del, dur, i0, amp_ge,amp_gi, pfreq, pdur,E_e, E_i

	USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
	  (nA) = (nanoamp)
  (mV) = (millivolt)
(umho) = (micromho)
}

PARAMETER {
	active = 0		
	active_pulses = 0	
	t_start = 0	(ms)	
	ksize = 20		
	imax = 2	(nA)

       del = 0	(ms)	
  dur = 10000	(ms)	
  
  E_e= 0 	  (mV)  
  E_i= -75 	  (mV)  

  amp_ge = 0	(umho)	
  amp_gi = 0	(umho)	
  pfreq = 10	(Hz)	
  pdur = 500	(ms)	
}

ASSIGNED {
	vcell	(mV)
	vr	(mV)
	vr10	(mV)
	vr50	(mV)
	ue	(mV)
	ue10	(mV)
	ue50	(mV)
	i_out	(nA)
	kernel[500]		
	previous_i[500]	(nA)	
	iz	(nA)
	ez	(mV)

       Npulses		
  pstart	(ms)	
  pstop		(ms)	
  g_e   (umho)          
  g_i   (umho)          
  g_eout (umho)         
  g_iout (umho)      	
}

INITIAL {
	FROM k = 0 TO 499 {
		previous_i[k]=0
	}
	vr = vcell
	vr10 = 10*vr
      vr50 = 50.4*vr
  g_e = 0
  g_i = 0
  Npulses = 0
  pstart = del
  pstop = pstart + pdur

	iz = 0
}

BREAKPOINT {
	SOLVE compensate
}


PROCEDURE compensate() {
	if (active && (i_out > imax || i_out < -imax)) {
		i_out = 0
	}

		vr = vcell
		if (active == 1 && t > t_start) {
			
			

			VERBATIM
			{
				int k,kmax;
				double *pi;
				double *kern;
				kmax = (int) ksize;
				pi = &(previous_i[kmax-1]);
				kern = &(kernel[kmax-1]);
				for(k = 1;k<kmax;k++) {
					vr-=(*pi=*(pi-1))*(*kern);
					pi--;
					kern--;
				}
				vr+=i_out*(*kern);
			}
			*previous_i = -i_out;
			ENDVERBATIM
		}
		vr10 = vr * 10
		vr50 = vr * 50.4
		ue = vcell-vr
		ue10 = ue * 10
	      ue50 = ue * 50

  if (active_pulses) {
if ((t >= pstart) && (t < dur)) { 



      g_e = amp_ge
      
      Npulses = Npulses + 1
      pstart = del + Npulses*1000/pfreq
      if (1000/pfreq <= pdur) { pstop = del + dur }

                               } 
    else if (t >= pstop) { 
      g_e = 0
      

      pstop = pstart + pdur
                         }
  }
  else { g_e = 0 
         
       }

iz = g_e * (vr - E_e) 

      g_eout = 50000.0 * g_e
      




}