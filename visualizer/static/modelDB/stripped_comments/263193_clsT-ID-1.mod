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
       SUFFIX clsT
       USEION na READ ena, nai, nao WRITE ina, nai, nao
       USEION k READ ek, ki, ko WRITE ik, ki, ko

       NONSPECIFIC_CURRENT il
       
       RANGE vLeftShift, AC, ACpotassium, timeLS, vLS0, vLS
       RANGE gnabar, gkbar, gna, gk, gnal, gkl, gl, el, INaKmax, ink, Kmko, Kmnai
       RANGE nai0, nao0, ki0, ko0
       RANGE qPump_0, qNa_0, qK_0, qGate_0 
       GLOBAL minf, hinf, ninf, mtau, htau, ntau, mLSinf, hLSinf, nLSinf, mLStau, hLStau, nLStau
  THREADSAFE 
}

PARAMETER {
        
        vLeftShift = 0.0 (mV) 
        timeLS = 0.0 (ms)
        vLS0 = 0.0 (mV)
          
        AC = 0.0  <0,1.0> 
        ACpotassium = 0.0  <0,1.0> 

        INaKmax = 9.09e-2       (mA/cm2) <0,1e6>




        ink0 = 0.010731         (mA/cm2) <0,1e6> 
        Kmnai =    10          (mM)    <0,1e6>
        Kmko =     3.5         (mM)    <0,1e6>



        
        nai0 = 20.0 (mM)
        nao0 =  154.0 (mM)
        ki0 =  150.0 (mM)
        ko0 =  6.0 (mM)


        gnabar = 0.12 (S/cm2) <0,1e9>
        gkbar = 0.036 (S/cm2) <0,1e9>
        gl = 0.0005 (S/cm2) <0,1e9>
        el = -59.9 (mV)
        gnal = 0.00025 (S/cm2) <0,1e9>
        gkl = 0.0001 (S/cm2) <0,1e9>

        
        qPump_0 = 1.9
        qNa_0 = 1.4
        qK_0 = 1.1
        qGate_0 = 3.0


}

STATE {

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
  LOCAL qPump, qNa, qK

  UNITSOFF
  qPump = qPump_0^((celsius - 20.0)/10.0)
  qNa = qNa_0^((celsius - 20.0)/10.0)
  qK = qK_0^((celsius - 20.0)/10.0)
  UNITSON

  SOLVE states METHOD cnexp 
  ink= qPump*INaKmax/(((1 + (Kmnai/nai))^3)*((1 + Kmko/ko)^2))

  gna = qNa*gnabar*(m*m*m*h*(1.0 - AC) + mLS*mLS*mLS*hLS*AC)
  ina = (gna + gnal)*(v - ena) + 3*ink
  gk = qK*gkbar*(n*n*n*n*(1.0 - ACpotassium) + nLS*nLS*nLS*nLS*ACpotassium )
  ik = (gk + gkl)*(v - ek) - 2*ink 
  il = gl*(v - el)
}




DERIVATIVE states {
        LOCAL nai_prime, ki_prime, qGate
        UNITSOFF
        qGate = qGate_0^((celsius - 20.0)/10.0)
        UNITSON

        rates(v)
        m' =  qGate*(minf-m)/mtau
        h' = qGate*(hinf-h)/htau
        n' = qGate*(ninf-n)/ntau

        mLS' =  qGate*(mLSinf-mLS)/mLStau
        hLS' = qGate*(hLSinf-hLS)/hLStau
        nLS' =  qGate*(nLSinf-nLS)/nLStau


        nai_prime = -(1e4)*4*ina/(FARADAY*diam)
        ki_prime = -(1e4)*4*ik/(FARADAY*diam)
        nai' = nai_prime
        nao' = -nai_prime
        ki' = ki_prime
        ko' = -ki_prime

}





PROCEDURE rates(v(mV)) {
      
      
      LOCAL  alpha, beta, sum
      

      UNITSOFF
      
      alpha = .1 * vtrap(-((v)+40),10)
      beta =  4 * exp(-((v)+65)/18)
      sum = alpha + beta
      mtau = 1/(sum)
      minf = alpha/sum
      
      alpha = .07 * exp(-((v)+65)/20)
      beta = 1 / (exp(-((v)+35)/10) + 1)
      sum = alpha + beta
      htau = 1/(sum)
      hinf = alpha/sum
      
      alpha = .01*vtrap(-(v+55),10) 
      beta = .125*exp(-(v+65)/80)
      sum = alpha + beta
      ntau = 1/(sum)
      ninf = alpha/sum

      if ( t < timeLS) {
               vLS = vLS0
       }else{
               vLS = vLeftShift
       }

      
      alpha = .1 * vtrap(-((v+vLS)+40),10)
      beta =  4 * exp(-((v+vLS)+65)/18)
      sum = alpha + beta
      mLStau = 1/(sum)
      mLSinf = alpha/sum
      
      alpha = .07 * exp(-((v+vLS)+65)/20)
      beta = 1 / (exp(-((v+vLS)+35)/10) + 1)
      sum = alpha + beta
      hLStau = 1/(sum)
      hLSinf = alpha/sum
      
      alpha = .01*vtrap(-((v+vLS)+55),10) 
      beta = .125*exp(-((v+vLS)+65)/80)
      sum = alpha + beta
      nLStau = 1/(sum)
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