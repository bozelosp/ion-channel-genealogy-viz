UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mol) = (1)
    (molar) = (mol/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
    (um) = (micrometer)
    FARADAY = (faraday) (coulombs)
}


NEURON {
       SUFFIX CLS
       USEION na READ ena, nai, nao WRITE ina, nai, nao
       USEION k READ ek, ki, ko WRITE ik, ki, ko

       NONSPECIFIC_CURRENT il
       
       RANGE vLS, vLS0, vLeftShift, AC, gnabar, gkbar, gna, gk, gnal, gkl, gl, el, ik, ina, INaKmax, ink, ink_last, Kmko, Kmnai, nai0, nao0, ki0, ko0, VolumeOut, VolumeIn, Area, time0, timeLS, ACpotassium
       GLOBAL minf, hinf, ninf, mtau, htau, ntau, mLSinf, hLSinf, nLSinf, mLStau, hLStau, nLStau
  THREADSAFE 
}

PARAMETER {
        CLSpotassium = 0 

        time0 = 0.0 (ms)
        timeLS = 200.0 (ms)
        vLS0 = 0.0 (mV)
        vLeftShift = 1.5 (mV)   
        AC = 1.0  <0,1.0> 
        ACpotassium = 0.0  <0,1.0> 

        INaKmax = 9.09e-2       (mA/cm2) <0,1e6>




        ink0 = 0.010731         (mA/cm2) <0,1e6> 
        Kmnai =    10          (mM)    <0,1e6>
        Kmko =     3.5         (mM)    <0,1e6>



        
        nai0 = 20.0 (mM)
        nao0 =  154.0 (mM)
        ki0 =  150.0 (mM)
        ko0 =  6.0 (mM)


        Area = 6.0 (um2)
        VolumeOut = 3.0 (um3)
        VolumeIn = 3.0 (um3)

        gnabar = 0.12 (S/cm2) <0,1e9>
        gkbar = 0.036 (S/cm2) <0,1e9>
        gl = 0.0005 (S/cm2) <0,1e9>
        el = -59.9 (mV)
        gnal = 0.00025 (S/cm2) <0,1e9>
        gkl = 0.0001 (S/cm2) <0,1e9>

}

STATE {
      time (ms)

      m
      h
      n
      mLS
      hLS
      nLS
      nai (mM)
      nao (mM)
      ki (mM)
      ko (mM)
}

ASSIGNED {
        v (mV)
        vLS (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

        ink (mA/cm2)
        ink_last (mA/cm2)

        gna (S/cm2)
        gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf
        hinf
        ninf
        mtau (ms) 
        htau (ms) 
        ntau (ms)

        mLSinf
        hLSinf
        nLSinf
        mLStau (ms) 
        hLStau (ms) 
        nLStau (ms) 

        diam (um)
        L (um)
}

INITIAL {
  v = -59.8137886
  vLS = vLS0
  time = time0  
  rates(v)
  m = minf
  h = hinf
  n = ninf
  mLS = mLSinf
  hLS = hLSinf
  nLS = nLSinf



  nai = nai0
  nao = nao0
  ki = ki0
  ko = ko0
  ink = ink0
}


BREAKPOINT {
  SOLVE states METHOD runge 
  
  ink= INaKmax/(((1 + (Kmnai/nai))^3)*((1 + Kmko/ko)^2))

  gna = gnabar*(m*m*m*h*(1.0 - AC) + mLS*mLS*mLS*hLS*AC)
  ina = (gna + gnal)*(v - ena) + 3*ink
  gk = gkbar*(n*n*n*n*(1.0 - ACpotassium) + nLS*nLS*nLS*nLS*ACpotassium )
  ik = (gk + gkl)*(v - ek) - 2*ink 
  il = gl*(v - el)
}




DERIVATIVE states {
        LOCAL nai_prime, ki_prime
        time' = 1.0

        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau

        mLS' =  (mLSinf-mLS)/mLStau
        hLS' = (hLSinf-hLS)/hLStau
        nLS' =  (nLSinf-nLS)/nLStau


        nai_prime = -(1e4)*4*ina/(FARADAY*diam)
        ki_prime = -(1e4)*4*ik/(FARADAY*diam)
        nai' = nai_prime
        nao' = -nai_prime
        ki' = ki_prime
        ko' = -ki_prime

}





PROCEDURE rates(v(mV)) {
      
      
      LOCAL  alpha, beta, sum, q10
      

      UNITSOFF
      q10 = 1.0 
      
      alpha = .1 * vtrap(-((v)+40),10)
      beta =  4 * exp(-((v)+65)/18)
      sum = alpha + beta
      mtau = 1/(q10*sum)
      minf = alpha/sum
      
      alpha = .07 * exp(-((v)+65)/20)
      beta = 1 / (exp(-((v)+35)/10) + 1)
      sum = alpha + beta
      htau = 1/(q10*sum)
      hinf = alpha/sum
      
      alpha = .01*vtrap(-(v+55),10) 
      beta = .125*exp(-(v+65)/80)
      sum = alpha + beta
      ntau = 1/(q10*sum)
      ninf = alpha/sum

      if ( time < timeLS) {
               vLS = vLS0
       }else{
               vLS = vLeftShift
       }

      
      alpha = .1 * vtrap(-((v+vLS)+40),10)
      beta =  4 * exp(-((v+vLS)+65)/18)
      sum = alpha + beta
      mLStau = 1/(q10*sum)
      mLSinf = alpha/sum
      
      alpha = .07 * exp(-((v+vLS)+65)/20)
      beta = 1 / (exp(-((v+vLS)+35)/10) + 1)
      sum = alpha + beta
      hLStau = 1/(q10*sum)
      hLSinf = alpha/sum
      
      alpha = .01*vtrap(-((v+vLS)+55),10) 
      beta = .125*exp(-((v+vLS)+65)/80)
      sum = alpha + beta
      nLStau = 1/(q10*sum)
      nLSinf = alpha/sum

}

FUNCTION vtrap(x,y) {  
       if (fabs(x/y) < 1e-6) {
               vtrap = y*(1 - x/y/2)
       }else{
               vtrap = x/(exp(x/y) - 1)
       }
}

UNITSON