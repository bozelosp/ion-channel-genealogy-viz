NEURON {
  POINT_PROCESS Fsquare

  RANGE del, dp, num, amp1, amp2
  POINTER x
}

PARAMETER {
  del = 0 (ms) <0, 1e9> 
  dp = 0 (ms) <0, 1e9> 
  num = 0 (1) 
  amp1 = 0 (1) 
  amp2 = 0 (1) 

}

ASSIGNED {
  x (1)
  on (1)
  tally (1) 
}

UNITSOFF
FUNCTION nonneg(x) {
  nonneg = x
  if (x<0) {
    nonneg = 0
  } else {
    nonneg = x
  }
}

INITIAL {
  on = 0
  x = 0
  
  del = nonneg(del)
  dp = nonneg(dp)
  num = nonneg(num)

  
  if (num*dp>0) {
    tally = num
    net_send(del,1) 
  }




}
UNITSON

NET_RECEIVE (w) {
  
  if (tally>0) { 
    if (flag == 1) { 
      if (on == 0) {
        on = 1 
      }
      x = amp1 
      
      net_send(dp, 2)
    }
    if (flag == 2) {
      x = amp2 
      tally = tally - 1
      if (tally>0) {
        net_send(dp, 1) 
      } else {
        net_send(dp, 3) 
      }
    }
  }
  if (flag == 3) { 
    on = 0 
    x = 0
  }
}