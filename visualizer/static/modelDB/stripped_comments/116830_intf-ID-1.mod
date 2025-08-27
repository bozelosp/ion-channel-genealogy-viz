NEURON {
  ARTIFICIAL_CELL INTF
  RANGE VAM, VNM, VGA, VGB, AHP            
  RANGE tauAM, tauNM, tauGA, tauGB, tauahp, ahpwt 
  RANGE VGBdel,tGB,VGBa,rebound,offsetGB   
  RANGE RMP,VTH,Vm,Vblock,refractory       
  RANGE taum,invl,oinvl,WINV,invlt         
  RANGE t0,tg,tGB,refrac                   
  RANGE nbur,tbur,cbur                     
  POINTER sop                              
  GLOBAL AMdec,NMdec,GAdec,GBdec           
  GLOBAL vdt,next,WEX,mg,RES,ESIN,Bb       
  GLOBAL tauGBGP,wGBGP,GPkd,Gn             
  GLOBAL EAM, ENM, EGA, EGB, spkht         
  GLOBAL prnum                             
}

PARAMETER {
  tauAM = 10 (ms)
  tauNM = 300 (ms)
  tauGA = 10 (ms)
  tauGB = 300 (ms)
  tauGBGP = 50 (ms) 
  taum =  10 (ms)
  invl =  100 (ms)
  WINV =  0
  wGBGP = 1 (ms) 
  GPkd  = 100    
  ahpwt = 0
  tauahp= 10 (ms)
  refrac = 5 (ms)
  VTH = -45      
  Vblock = -20   
  vdt = 0.1      
  mg = 1         
  sop=0
  AMdec=1       
  NMdec=1
  GAdec=1
  GBdec=1
  nbur=1
  tbur=2
  VGBdel=30
  rebound=0.01 
  offsetGB=0
  RMP=-65
  EAM = 65
  ENM = 90
  EGA = -15
  EGB = -30
  spkht = 50
  prnum = -1
}

ASSIGNED {
  VAM
  VNM
  VGA
  VGB
  VGBa
  AHP
  t0(ms)
  tGB(ms)
  tg(ms)
  refractory
  next
  WEX
  RES
  ESIN
  Gn
  Bb
  cbur
  Vm
  invlt
  oinvl
}



CONSTRUCTOR {
  VERBATIM 
  { int lid,lty;
    if (ifarg(2)) { lid=(int) *getarg(2); } else { lid= UINT_MAX; }
    if (ifarg(3)) { lty=(int) *getarg(3); } else { lty= -1; }
    _p_sop = (double*)ecalloc(1, sizeof(id0));
    ip = IDP;
    ip->id=lid; ip->type=lty; 
    ip->invl0 = ip->record = ip->jitter = ip->input = 0; 
    ip->rve=-1;
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
  free(IDP);
  }
  ENDVERBATIM
}

INITIAL {
  VAM = 0
  VNM = 0
  VGA = 0
  VGB = 0
  VGBa = 0
  t0 = t
  tGB = t
  tg = 0
  offsetGB=0
  AHP=0
  invlt = -1
  VERBATIM
  jtpt=0;    
  errflag=0;
  ENDVERBATIM
  refractory = 0 
  
  if (vinflag()) { randspk() net_send(next,2) }
  if (recflag()) { recini() } 
}


NET_RECEIVE (wAM,wNM,wGA,wGB,wflg) { LOCAL tmp,rflg,wrec,id,jflg,iflg,invlflg,ty
  INITIAL { wNM=wNM wGA=wGA wGB=wGB wflg=0}
  
  VERBATIM
  ip = IDP;  
  _lrflg=(double)ip->record; _ljflg=(double)ip->jitter; _liflg=(double)ip->input; 
  _linvlflg=(double)ip->invl0; _lwrec=(double)ip->wrec; _lid=(double)ip->id; 
  _lty=(double)ip->type; 
  ENDVERBATIM
  if (flag==4) { 
    cbur=cbur-1  
    if (cbur>0) { net_send(tbur,4) 
    } else { net_send(refrac-AHP/10, 3) }
    if (jflg) { tmp= t+jitter()/10 } else { tmp=t }
    net_event(tmp)
    if (rflg) { recspk(tmp) }
    if (wrec) { wrecord(tmp,-1,0) }
  
  
  
  } else if (flag==0 && wGB==0 && wflg==1) {
    VERBATIM
    ip->input=1;
    ENDVERBATIM
    iflg=1 
    wflg=2 
    randspk() 
    net_send(next,2)
  } else if (flag==0 && wGB==0 && wflg==2) { 
    VERBATIM
    ip->input=0;
    ENDVERBATIM
    iflg=0 
    wflg=1  
  } else {  
    if (rflg) { record() }
    
    if (VAM>hoc_epsilon)  { VAM = VAM*EXP(-(t - t0)/tauAM) } else { VAM=0 } 
    if (VNM>hoc_epsilon)  { VNM = VNM*EXP(-(t - t0)/tauNM) } else { VNM=0 } 
    if (VGA< -hoc_epsilon){ VGA = VGA*EXP(-(t - t0)/tauGA) } else { VGA=0 } 
    VGB = esinr(t-tGB) 
    if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } 
    
    Vm = VAM+VNM+VGA+VGB+AHP 
    if (Vm>100||Vm<-60){ pid() 
      printf("WARN
    if (flag==0) { 
      
      if (wAM>0 && VAM<EAM) { 
        tmp = wAM*(1-Vm/EAM)
        VAM = VAM + tmp
        if (wrec) { wrecord(t,0,tmp) }
      }
      
      if (wNM>0 && VNM<ENM) { rates(RMP+Vm)
        tmp = wNM*Bb*(1-Vm/ENM) 
        VNM = VNM + tmp
        if (wrec) { wrecord(t,1,tmp) }
      } 
      if (VNM>1.2*ENM) { pid() 
        
        printf("**** ERR
      }
      
      if (wGA>0 && VGA>EGA) { 
        tmp = wGA*(1-Vm/EGA) 
        VGA = VGA - tmp
        if (wrec) { wrecord(t,2,tmp) }
      }
      if (wGB>0) { net_send(VGBdel,5) } 
      if (invlflg) { 
        Vm = RMP+VAM+VNM+VGA+VGB+AHP
        if (invlt==-1) { 
          if (Vm>RMP) {
            oinvl=invl
            invlt=t
            net_send(invl,1) 
          }
        } else {
          tmp=shift(Vm)
          if (tmp!=0)  {
            net_move(tmp) 
            if (id()<prnum) {
               pid() printf("**** MOVE t=%g to %g Vm=%g %g,%g\n",t,tmp,Vm,invlt,oinvl) }
          }
        }      
      }
    } else if (flag==5) { 
      offsetGB = VGB 
      
      
      wflg=wflg*EXP(-(t-tGB)/tauGBGP)+wGBGP 
      coop(wflg)               
      
      if (VGB>EGB) { 
        tmp = wGB*(1-Vm/EGB)*Gn 
        VGB = VGB - tmp
        if (wrec) { wrecord(t,3,tmp) }
      }
      VGBa= VGB
      tGB=t 
    
    } else if (flag==2) { 
      if (iflg==0) { flag=-1 
      } else { 
        randspk() 
        net_send(next,2) 
        if (WEX<0) { 
          net_event(t)   
          if (rflg) { recspk(t) }
          if (wrec) { wrecord(t,-1,0) }
        } else {
          tmp = WEX*(1-Vm/EAM)
          VAM = VAM + tmp
          if (wrec) { wrecord(t,0,tmp) }
        }
      }
    } else if (flag==1) { 
      
      if (WINV<0) { 
        net_event(t)   
        if (rflg) { recspk(t) }
        if (wrec) { wrecord(t,-1,0) }
      } else {
        tmp = WINV*(1-Vm/EAM)
        VAM = VAM + tmp 
        if (wrec) { wrecord(t,0,tmp) }
      }
      oinvl=invl
      invlt=t
      net_send(invl,1) 
    } else if (flag==3) { 
      refractory = 0 
    }
    
    Vm = VAM+VNM+VGA+VGB+RMP+AHP
    if (refractory==0 && (Vm>VTH && Vm<Vblock)) {
      AHP = AHP - ahpwt
      if (jflg) { tmp=t+jitter() } else { tmp=t }
      net_event(tmp)
      if (rflg) { recspk(tmp) }
      if (wrec) { wrecord(tmp,-2,0) }
      VAM=VAM*AMdec VNM=VNM*NMdec
      VGA=VGA*GAdec VGB=VGB*GBdec
      if (nbur>1) { cbur=nbur-1 net_send(tbur,4) 
      } else { net_send(refrac-AHP/10, 3) } 
      refractory = 1
    }
    t0 = t
  }
}



PROCEDURE randspk () {
  VERBATIM 
  ip=IDP;  
  if (ip->rvi > ip->rve) { ip->input=0;
  } else { 
    WEX=wsp[ip->rvi];
    next=vsp[ip->rvi++]-t; 
  }
  ENDVERBATIM
  
}


VERBATIM
void val (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  vii[4]=esinr(xx-tGB);
  vii[5]=AHP*EXP(-(xx - ta)/tauahp);
  vii[6]=vii[1]+vii[2]+vii[3]+vii[4]+vii[5];
}
ENDVERBATIM


PROCEDURE record () {
  VERBATIM {
  int k; double ti;
  vp = SOP;
  if (tg>=t) return 0;
  if (vp->p >= vp->size) { if (errflag) return 0; 
    printf("**** WARNING out of recording room for INTF type%d id%d at %g****\n",IDP->type,IDP->id,t);
    printf("**************** WARNING: No further WARNINGS ****************\n");
    errflag=1; return 0; }
  for (ti=tg;ti<=t && vp->p < vp->size;ti+=vdt,vp->p++) { 
    val(ti,tg);  
    vp->vvo[0][vp->p]=ti;
    for (k=1;k<NSV;k++) if (vp->vvo[k]!=0) { 
      vp->vvo[k][vp->p]=vii[k]+RMP;
    }
  }
  tg=t;
  }
  ENDVERBATIM
}


PROCEDURE recspk (x) {
  VERBATIM { int k;
  vp = SOP;
  record();
  if (vp->p >= vp->size || vp->vvo[6]==0) return 0; 
  vp->vvo[0][vp->p-1]=_lx;
  vp->vvo[6][vp->p-1]=spkht; 
  tg=_lx;
  }
  ENDVERBATIM
}


PROCEDURE recclr () {
  VERBATIM 
  {int k;
  if (IDP->record) {
    if (SOP!=nil) {
      vp = SOP;
      vp->size=0; vp->p=0;
      for (k=0;k<NSV;k++) { vp->vv[k]=nil; vp->vvo[k]=nil; }
    } else printf("INTF recclr ERR: nil pointer\n");
  }
  IDP->record=0;
  }
  ENDVERBATIM 
}


PROCEDURE recfree () {
  VERBATIM
  if (SOP!=nil) {
    free(SOP);
    SOP=nil;
  } else printf("INTF recfree ERR: nil pointer\n");
  IDP->record=0;
  ENDVERBATIM
}





PROCEDURE initvspks () {
  VERBATIM 
  {int max, i, err=0; 
    double *iv, last;
    if (! ifarg(1)) {printf("Return initvspks(ivspks,vspks,wvspks)\n"); return 0.;}
    ip=IDP;
    i = vector_arg_px(1, &iv);
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("initvspks ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("initvspks ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("initvspks ERR: 3rd vec is of different size\n");}
    ip->vinflg=1;
    for (i=0; i<max && (int)iv[i] != ip->id ; i++); 
    if (i==max) { 
      printf("initvspks WARN: %d not found in ivspks\n",ip->id); 
      ip->vinflg=0; ip->rve=-1;
      return(0.); 
    }
    ip->rvb=ip->rvi=i;
    last=vsp[i++];
    for (; i<max && (int)iv[i] == ip->id ; i++) { 
      if (vsp[i]<=last) { err=1; 
        printf("initvspks ERR: nonmonotonic for cell#%d: %g %g\n",ip->id,last,vsp[i]); }
      last=vsp[i];
    }
    ip->rve=i-1;
    if (err) { ip->rve=0; hoc_execerror("",0); }
  }
  ENDVERBATIM
}



PROCEDURE initjitter () {
  VERBATIM 
  {int max, i, err=0;
    jtpt=0;
    if (! ifarg(1)) {printf("Return initjitter(vec)\n"); return(0.);}
    max=vector_arg_px(1, &jsp);
    if (max==0) {err=1; printf("initjitter ERR: vec not initialized\n");}
    for (i=0; i<max; i++) if (jsp[i]<=0) {err=1;
      printf("initjitter ERR: vec should be >0: %g\n",jsp[i]);}
    if (err) { jsp=nil; jitmax=0.; return(0.); }
    if (max != jitmax) {
      printf("WARNING: resetting jitmax_INTF to %d\n",max); jitmax=max; }
  }
  ENDVERBATIM
}



PROCEDURE initinvl () {
  VERBATIM 
  printf("initinvl() NOT BEING USED\n"); return(0.);
  ENDVERBATIM
}


FUNCTION invlflag () {
  VERBATIM
  ip=IDP;
  if (ip->invl0==1 && invlp==nil) { 
    printf("INTF invlflag ERR: pointer not initialized\n"); hoc_execerror("",0); 
  }
  _linvlflag= (double)ip->invl0;
  ENDVERBATIM
}


FUNCTION shift (vl) { 
  VERBATIM   
  double expand, tmp, min, max;

  if ((t<(invlt-invl)+invl/2) && invlt != -1) { 
    _lshift=0.;  
  } else {
    expand = -(_lvl-(-65))/20; 
    if (expand>1.) expand=1.; if (expand<-1.) expand=-1.;
    if (expand>0.) { 
      max=1.5*invl;
      tmp=oinvl+0.8*expand*(max-oinvl); 
    } else {
      min=0.5*invl; 
      tmp=oinvl+0.8*expand*(oinvl-min); 
    }
    if (invlt+tmp<t+2) { 
      _lshift=0.;
    } else {
      oinvl=tmp; 
      _lshift=invlt+oinvl;
    }
  }
  ENDVERBATIM
}


PROCEDURE recini () {
  VERBATIM 
  { int k;
  if (SOP==nil) { 
    printf("INTF record ERR: but pointer not initialized\n"); hoc_execerror("",0); 
  } else {
    vp = SOP;
    vp->p=0;
    
    for (k=0;k<NSV;k++) if (vp->vvo[k]!=0) vector_resize(vp->vv[k], vp->size);
  }}
  ENDVERBATIM
}


PROCEDURE fini () {
  VERBATIM 
  {int k;
  
  IDP->rvi=IDP->rvb;  
  if (IDP->record) {
    record(); 
    for (k=0;k<NSV;k++) if (vp->vvo[k]!=0) { 
      vector_resize(vp->vv[k], vp->p);
    }
  }}
  ENDVERBATIM
}



PROCEDURE chk () {
  VERBATIM 
  {int i,lfg;
  ip=IDP;
  printf("ID:%d; typ: %d; rec:%d wrec:%d inp:%d jit:%d invl:%d\n",ip->id,ip->type,ip->record,ip->wrec,ip->input,ip->jitter,ip->invl0);
  if (ifarg(1)) lfg=(int) *getarg(1); else lfg=0;
  if (lfg==1) {
    if (SOP!=nil) {
      vp = SOP;
      printf("p %d size %d tg %g\n",vp->p,vp->size,tg);
      for (i=0;i<NSV;i++) printf("%d %p %p;",i,vp->vv[i],vp->vvo[i]);
    } else printf("Recording pointers not initialized");
  }
  if (lfg==2) { 
    printf("Global vectors for input and jitter: \n");
    if (vsp!=nil) printf("VSP: %p (%d/%d-%d)\n",vsp,ip->rvi,ip->rvb,ip->rve); else printf("no VSP\n");
    if (jsp!=nil) printf("JSP: %p (%d/%d)\n",jsp,jtpt,jitmax); else printf("no JSP\n");
  }
  if (lfg==3) { 
    if (vsp!=nil) { printf("VSP: (%d/%d-%d)\n",ip->rvi,ip->rvb,ip->rve); 
      for (i=ip->rvb;i<=ip->rve;i++) printf("%d:%g  ",i,vsp[i]);
      printf("\n");
    } else printf("no VSP\n");
  }
  if (lfg==4) {  
  }
  if (lfg==5) { 
    printf("wwpt %d wwsz %d\n WW vecs: ",wwpt,wwsz);
    for (i=0;i<NSW;i++) printf("%d %p %p;",i,ww[i],wwo[i]);
  }}
  ENDVERBATIM
}


FUNCTION pid () {
  VERBATIM 
  printf("INTF%d(%d) ",IDP->id,IDP->type);
  _lpid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION id () {
  VERBATIM 
  _lid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION type () {
  VERBATIM 
  _ltype = (double)IDP->type;
  ENDVERBATIM
}


PROCEDURE initrec () {
  VERBATIM 
  {int i; void *vv;
  name = gargstr(1);
  if (SOP==nil) { 
    IDP->record=1;
    SOP = (vpt*)ecalloc(1, sizeof(vpt));
    SOP->size=0;
  }
  if (IDP->record==0) {
    recini();
    IDP->record=1;
  }
  vp = SOP;
  i=(int)varnum();
  if (i==-1) {printf("INTF record ERR %s not recognized\n",name); hoc_execerror("",0); }
  vp->vv[i]=vector_arg(2);
  vector_arg_px(2, &(vp->vvo[i]));
  if (vp->size==0) { vp->size=(unsigned int)vector_buffer_size(vp->vv[i]);
  } else if (vp->size != (unsigned int)vector_buffer_size(vp->vv[i])) {
    printf("INTF initrec ERR vectors not all same size: %d vs %d",vp->size,vector_buffer_size(vp->vv[i]));
    hoc_execerror("", 0); 
  }} 
  ENDVERBATIM
}



FUNCTION varnum () { LOCAL i
  i=-1
  VERBATIM
  if (strcmp(name,"time")==0)      { _li=0.;
  } else if (strcmp(name,"VAM")==0) { _li=1.;
  } else if (strcmp(name,"VNM")==0) { _li=2.;
  } else if (strcmp(name,"VGA")==0) { _li=3.;
  } else if (strcmp(name,"VGB")==0) { _li=4.;
  } else if (strcmp(name,"AHP")==0) { _li=5.;
  } else if (strcmp(name,"V")==0) { _li=6.;
  } else if (strcmp(name,"VM")==0) { _li=6.; 
  }
  ENDVERBATIM
  varnum=i
}


PROCEDURE vecname () {
  VERBATIM
  int i; 
  i = (int)*getarg(1);
  if (i==0)      printf("time\n");
  else if (i==1) printf("VAM\n");
  else if (i==2) printf("VNM\n");
  else if (i==3) printf("VGA\n");
  else if (i==4) printf("VGB\n");
  else if (i==5) printf("AHP\n");
  else if (i==6) printf("V\n");
  ENDVERBATIM
}


PROCEDURE rebuild () { LOCAL s0,w0,wwaz,ii,wflg,tmp
  VERBATIM
  int ii,jj; double i0,tstop;
  ip = IDP; vp=SOP; 
  ip->record=1; ip->input=0; i0=wwo[1][0];
  initmodel(); _lwwaz=wwaz; _lii=0.;
  tstop=(ifarg(1))?(*getarg(1)):0.;
  ENDVERBATIM
  WHILE (ii<wwaz) { 
    VERBATIM
    int ii=(int)_lii;
    if (wwo[1][ii]!=i0) {
      printf("ERROR wrong id at %d %g not %g\n",ii,wwo[1][ii],i0); hoc_execerror("", 0);}
    t=wwo[0][ii]; _ls0=wwo[2][ii]; _lw0=wwo[3][ii];
    ENDVERBATIM
    record()
    
    if (VAM>hoc_epsilon)  { VAM = VAM*EXP(-(t - t0)/tauAM) } else { VAM=0 } 
    if (VNM>hoc_epsilon)  { VNM = VNM*EXP(-(t - t0)/tauNM) } else { VNM=0 } 
    if (VGA< -hoc_epsilon){ VGA = VGA*EXP(-(t - t0)/tauGA) } else { VGA=0 } 
    VGB = esinr(t-tGB) 
    if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } 
    if (s0==0) { VAM = VAM + w0 }
    if (s0==1) { VNM = VNM + w0 }
    if (s0==2) { VGA = VGA - w0 }
    if (s0==3) { 
      offsetGB = VGB 
      VGB = VGB - w0
      VGBa= VGB
      tGB=t
    }
    if (s0==-2) { 
      recspk(t) 
      AHP = AHP - ahpwt
      VAM=VAM*AMdec VNM=VNM*NMdec
      VGA=VGA*GAdec VGB=VGB*GBdec
    }
    if (s0==-1) { recspk(t) }
    t0 = t
    ii = ii+1
  }
  if (tstop>0) {
    VERBATIM
    t=tstop;  
    ENDVERBATIM
    record()
  }
}



PROCEDURE initwrec () {
  VERBATIM 
  {int k;
  if (! ifarg(4)) { 
    wwsz=WSZ;
    wf1 = hoc_obj_file_arg(1);
    wf2 = hoc_obj_file_arg(2);
  } else if (! ifarg(8)) { 
    if (NSW!=4) { 
      printf("INTF initwrec ERR w-vecs compiled for 4 args\n");
      hoc_execerror("",0); }
    IDP->wrec=1;
    for (k=0;k<NSW;k++) {
      ww[k]=vector_arg(k+1);
      wwaz=vector_arg_px(k+1, &(wwo[k]));
    }
    if (wwsz==0) wwsz=(unsigned int)vector_buffer_size(ww[0]); 
    for (k=0;k<NSW;k++) if (wwsz!=(unsigned int)vector_buffer_size(ww[k])) {
      printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,vector_buffer_size(ww[k]));
    }
  } else { 
    if (NSW!=8) { 
      printf("INTF initwrec ERR w-vecs compiled for 8 args\n");
      hoc_execerror("",0); }
    IDP->wrec=1;
    for (k=0;k<NSW;k++) {
      ww[k]=vector_arg(k+1);
      wwaz=vector_arg_px(k+1, &(wwo[k]));
    }
    if (wwsz==0) wwsz=(unsigned int)vector_buffer_size(ww[0]); 
    for (k=0;k<NSW;k++) if (wwsz!=(unsigned int)vector_buffer_size(ww[k])) {
      printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,vector_buffer_size(ww[k]));
    }
  }}
  ENDVERBATIM
}


PROCEDURE wrecord (t,s0,w0) {
  VERBATIM {
  int k; double id = (double)IDP->id;
  if (wwpt >= wwsz) { 
    wwpt=0;
    fprintf(wf1,"
    fwrite(&wwt,sizeof(float),WSZ,wf2);  
    fwrite(&wwi,sizeof(int),WSZ,wf2);  
    fwrite(&wws,sizeof(char),WSZ,wf2);  
    fwrite(&www,sizeof(float),WSZ,wf2);  
  } 
  
  wwt[wwpt]=(float)_lt; 
  wwi[wwpt]=(unsigned int)IDP->id; 
  wws[wwpt]=(char)_ls0; 
  www[wwpt]=(float)_lw0; 
  wwpt++;
  }
  ENDVERBATIM
}

FUNCTION wrec () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->wrec = (short) *getarg(1);
  _lwrec=(double)ip->wrec;
  ENDVERBATIM
}

FUNCTION wwszset () {
  VERBATIM
  if (ifarg(1)) wwsz = (short) *getarg(1);
  _lwwszset=(double)wwsz;
  ENDVERBATIM
}


FUNCTION wwfree () {
  VERBATIM
  int k;
  IDP->wrec=0;
  wwsz=0; wwpt=0;
  for (k=0;k<NSW;k++) { ww[k]=nil; wwo[k]=nil; }
  ENDVERBATIM
}



FUNCTION jitter () {
  if (jitmax>0 && jtpt>=jitmax) {  jtpt=0
    printf("Warning, cycling through jitter vector at t=%g\n",t) }
  if (jitmax>0) {
    VERBATIM 
    _ljitter = jsp[jtpt++];
    ENDVERBATIM
  } else { jitter=0 }
}


PROCEDURE global_init () {
  VERBATIM
  int k;
  if (wwo[0]!=0) { 
    printf("Initializing global ww vectors\n");
    wwpt=0;
    for (k=0;k<NSW;k++) vector_resize(ww[k], wwsz);
  }
  ENDVERBATIM
}

PROCEDURE global_fini () {
  VERBATIM
  {int k;
  if (IDP->wrec) {
    if (wwo[0]!=0) { 
      for (k=0;k<NSW;k++) vector_resize(ww[k], wwpt);
    } else {
      fprintf(wf1,"
      fwrite(&wwt,sizeof(float),wwpt,wf2);  
      fwrite(&wwi,sizeof(int),wwpt,wf2);  
      fwrite(&wws,sizeof(char),wwpt,wf2);  
      fwrite(&www,sizeof(float),wwpt,wf2);  
      printf("Closing file with wwpt=%d at location %ld\n",wwpt,ftell(wf2));
      fclose(wf1); fclose(wf2);
    }
  } else {
  printf("WARNING: global_fini() called from %d:%d with no wrec pointers\n",IDP->type,IDP->id);
  }}
  ENDVERBATIM
}


FUNCTION fflag () { fflag=1 }
FUNCTION thresh () { thresh=VTH-RMP }


FUNCTION recflag () { 
  VERBATIM
  _lrecflag= (double)IDP->record;
  ENDVERBATIM
}


FUNCTION vinflag () {
  VERBATIM
  ip=IDP;
  if (ip->vinflg==0 && vsp==nil) { 
  } else if (ip->vinflg==1 && ip->rve==-1) {
    printf("INTF vinflag ERR: pointer not initialized\n"); hoc_execerror("",0); 
  } else if (ip->rve >= 0) { 
    if (vsp==nil) {
      printf("INTF vinflag ERR1: pointer not initialized\n"); hoc_execerror("",0); 
    }
    ip->rvi=ip->rvb;
    ip->input=1;
  }
  _lvinflag= (double)ip->vinflg;
  ENDVERBATIM
}


FUNCTION jitset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->jitter = (short) *getarg(1);
  _ljitset=(double)ip->jitter;
  ENDVERBATIM
}


FUNCTION invlset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->invl0 = (short) *getarg(1);
  _linvlset=(double)ip->invl0;
  ENDVERBATIM
}


FUNCTION vinset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->vinflg = (short) *getarg(1);
  if (ip->vinflg==1) {
    ip->input=1;
    ip->rvi = ip->rvb;
  }
  _lvinset=(double)ip->vinflg;
  ENDVERBATIM
}


PROCEDURE EXPo (x) {
  TABLE RES FROM -20 TO 0 WITH 5000
  RES = exp(x)
}

FUNCTION EXP (x) {
  EXPo(x)
  EXP = RES
}

FUNCTION esinr (x) {
  ESINo(PI*x/tauGB)
  if        (x<tauGB)   {    esinr= (VGBa-offsetGB)*ESIN +offsetGB
  } else if (x>2*tauGB) {    esinr= 0
  } else {                   esinr= rebound*VGBa*ESIN }
}

PROCEDURE ESINo (x) {
  TABLE ESIN FROM 0 TO 2*PI WITH 3000 
  ESIN = sin(x)
}

PROCEDURE rates(vv) {
  TABLE Bb DEPEND mg FROM -100 TO 50 WITH 300
  
  Bb = 1 / (1 + exp(0.062 (/mV) * -vv) * (mg / 3.57 (mM)))
}

PROCEDURE coop (x) {
  TABLE Gn DEPEND GPkd FROM 0 TO 10 WITH 100
  
  Gn = (x^4)/(x^4+GPkd) 
}