NEURON {
  ARTIFICIAL_CELL INTF
  RANGE VAM, VNM, VGA, VGB, AHP            
  RANGE tauAM, tauNM, tauGA, tauGB, tauahp, ahpwt 
  RANGE VGBdel,tGB,VGBa,rebound,rebob,offsetGB   
  RANGE RMP,VTH,Vm,Vblock,refractory       
  RANGE taum,invl,oinvl,WINV,invlt         
  RANGE t0,tg,twg,tGB,refrac,Vbrefrac      
  RANGE nbur,tbur,cbur,AHP2REF,WEX         
  POINTER sop                              
  GLOBAL AMdec,NMdec,GAdec,GBdec           
  GLOBAL vdt,next,mg,RES,ESIN,Bb,Psk   
  GLOBAL tauGBGP,wGBGP,GPkd,Gn             
  GLOBAL EAM, ENM, EGA, EGB, spkht         
  GLOBAL prnum,wwwid,wwht,nsw,rebeg        
  GLOBAL subsvint
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
  Vbrefrac = 20 (ms)
  wwwid = 10
  wwht = 10
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
  VGBdel=0
  rebound=0.01 
  offsetGB=0
  RMP=-65
  EAM = 65
  ENM = 90
  EGA = -15
  EGB = -30
  spkht = 50
  prnum = -1
  nsw=0
  AHP2REF=0.1
  rebeg=0
  subsvint=0
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
  twg(ms)
  refractory
  next
  WEX
  RES
  ESIN
  Gn
  Bb
  Psk
  cbur
  Vm
  invlt
  oinvl
  rebob
}



CONSTRUCTOR {
  VERBATIM 
  { int lid,lty,lco;
    if (ifarg(1)) { lid=(int) *getarg(1); } else { lid= UINT_MAX; }
    if (ifarg(2)) { lty=(int) *getarg(2); } else { lty= -1; }
    if (ifarg(3)) { lco=(int) *getarg(3); } else { lco= -1; }
    _p_sop = (double*)ecalloc(1, sizeof(id0));
    ip = IDP;
    ip->id=lid; ip->type=lty; ip->col=lco; 
    ip->invl0 = ip->record = ip->jitter = ip->input = 0; 
    ip->vbr=0;
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
  twg = 0
  offsetGB=0
  AHP=0
  rebob=-1e9
  invlt = -1
  VERBATIM
  jtpt=0;    
  errflag=0;
  ENDVERBATIM
  refractory = 0 
  
  if (vinflag()) { randspk() net_send(next,2)}
  if (recflag()) { recini() } 
  rebeg=0 
}


NET_RECEIVE (wAM,wNM,wGA,wGB,wflg) { LOCAL tmp
 INITIAL { wNM=wNM wGA=wGA wGB=wGB wflg=0}
  
VERBATIM
  ip = IDP;  

ENDVERBATIM
VERBATIM
  if (ip->dbx>2) 
ENDVERBATIM
{ 
    pid() 
    printf("DB0
    if (flag==0) { printf(" (%g %g %g %g %g)",wAM,wNM,wGA,wGB,wflg) }
    printf("\n")
  }
  
  if (flag==4) { 
    cbur=cbur-1  
    if (cbur>0) { 
      net_send(tbur,4) 
    } else { 
      net_send(refrac-AHP*AHP2REF, 3) 
    }
    tmp=t
VERBATIM
    if (ip->jitter) 
ENDVERBATIM
{ tmp= t+jitter()/10 } 
    net_event(tmp)
VERBATIM
    if (ip->dbx>0) 
ENDVERBATIM
{ pid() printf("DBA
VERBATIM
    if (ip->record) 
ENDVERBATIM
{ recspk(tmp) } 
VERBATIM
    if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
    
  
  
  
  
  } else if (flag==0 && wGB==0 && wflg==1) {
VERBATIM
    ip->input=1; 

ENDVERBATIM
    wflg=2 
    randspk() 
    net_send(next,2)
  } else if (flag==0 && wGB==0 && wflg==2) { 
VERBATIM
    ip->input=0; 

ENDVERBATIM
    wflg=1  
  } else {  
VERBATIM
    if (ip->record) 
ENDVERBATIM
{ record() } 
VERBATIM
    if (ip->wrec) 
ENDVERBATIM
{ wrecord(1e9) } 
    
    if (VAM>hoc_epsilon)  { VAM = VAM*EXP(-(t - t0)/tauAM) } else { VAM=0 } 
    if (VNM>hoc_epsilon)  { VNM = VNM*EXP(-(t - t0)/tauNM) } else { VNM=0 } 
    if (VGA< -hoc_epsilon){ VGA = VGA*EXP(-(t - t0)/tauGA) } else { VGA=0 } 
    if (VGBdel>0) {
      VGB = esinr(t-tGB) 
    } else {
      if (VGB< -hoc_epsilon){ 
        VGB = VGB*EXP(-(t - t0)/tauGB) } else { VGB=0 }
    }      
    if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } 
    t0 = t 
    
    Vm = VAM+VNM+VGA+VGB+AHP 
    if (Vm>100||Vm<-60){ pid() 
      printf("ERR_A
      stoprun=1
    }
    if (flag==0) { 
      
      if (wAM>0) {
        if (rebob!=1e9 && rebob!=-1e9) {
VERBATIM
          cbur=floor(rebound*rebob/EGB); 

ENDVERBATIM
VERBATIM
          if (ip->dbx==-1) 
ENDVERBATIM
{ pid() printf("C
          net_send(tbur,4) 
          rebob=1e9
        }
        if (VAM<EAM) {
          tmp = wAM*(1-Vm/EAM)
          VAM = VAM + tmp
        }
      }
      
      if (wNM>0 && VNM<ENM) { rates(RMP+Vm)
        tmp = wNM*Bb*(1-Vm/ENM) 
        VNM = VNM + tmp
        
      } 
      if (VNM>1.2*ENM) { pid() 
        
        printf("ERR_B
        stoprun=1
      }
      
      if (wGA>0 && VGA>EGA) { 
        tmp = wGA*(1-Vm/EGA) 
        VGA = VGA - tmp
        
      }
      if (wGB>1e-6) {
        if (VGBdel>0) { net_send(VGBdel,5)  
        } else { 
          
          wflg=wflg*EXP(-(t-tGB)/tauGBGP)+wGBGP 
          coop(wflg)               
          tmp = wGB*(1-Vm/EGB)*Gn
          VGB = VGB - tmp
          if (VGB<rebob && rebob!=1e9 && rebob!=-1e9) { rebob=VGB }
        }
      }
VERBATIM
      if (ip->invl0) 
ENDVERBATIM
{ 
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
        
      }
      VGBa= VGB
      tGB=t 
    
    } else if (flag==2) { 
VERBATIM
      if (ip->input==0) 
ENDVERBATIM
{ flag=-1 } 
      if (flag==2) { 
{pid() printf("DBBa
        if (WEX<0) { 
          net_event(t)   
VERBATIM
          if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBB
          if (WEX<-1) { cbur=-WEX  net_send(tbur,4) }
VERBATIM
          if (ip->record) 
ENDVERBATIM
{ recspk(t) } 
VERBATIM
          if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
          
        } else if (WEX>0) {
          tmp = WEX*(1-Vm/EAM)
          VAM = VAM + tmp
          
        }
        randspk() 
        if (next>0) { net_send(next,2) }
      }
    } else if (flag==1) { 
      
      if (WINV<0) { 
        net_event(t)   
VERBATIM
        if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBC
VERBATIM
        if (ip->record) 
ENDVERBATIM
{ recspk(t) } 
VERBATIM
        if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
        
      } else {
        tmp = WINV*(1-Vm/EAM)
        VAM = VAM + tmp 
        
      }
      oinvl=invl
      invlt=t
      net_send(invl,1) 
    } else if (flag==3) { 
      refractory = 0 
    }
    Vm = VAM+VNM+VGA+VGB+RMP+AHP
    if (refractory==0 && Vm>VTH) {
VERBATIM
      if (!ip->vbr && Vm>Vblock) return; 

ENDVERBATIM
      AHP = AHP - ahpwt
      tmp=t
VERBATIM
      if (ip->jitter) 
ENDVERBATIM
{ tmp= t+jitter() }  
      net_event(tmp)
VERBATIM
      if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBD
VERBATIM
      if (ip->record) 
ENDVERBATIM
{ recspk(tmp) } 
VERBATIM
      if (ip->wrec) 
ENDVERBATIM
{ wrecord(tmp) } 
      
      VAM=VAM*AMdec VNM=VNM*NMdec
      VGA=VGA*GAdec VGB=VGB*GBdec
      if (nbur>1) { 
        cbur=nbur-1 net_send(tbur,4) 
VERBATIM
        return; 

ENDVERBATIM
      } else if (rebob==1e9) { rebob=0 }
      refractory = 1
VERBATIM
      if (ip->vbr && Vm>Vblock) 
ENDVERBATIM
{ 
        net_send(Vbrefrac,3) 
VERBATIM
        if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBE
VERBATIM
        return; 

ENDVERBATIM
      }
      net_send(refrac-AHP*AHP2REF, 3) 
    }
  }
}



PROCEDURE randspk () {
  VERBATIM 
  ip=IDP;  
  if (ip->rvi > ip->rve) { 
    ip->input=0;           
    next=-1.;
  } else { 
    
    while ((next=vsp[ip->rvi++]-t)<=1e-6) if (ip->rvi > ip->rve) { 
      printf("randspk() ERRA: "); chk(2.); hxe(); }
    WEX=wsp[ip->rvi-1]; 
    if (ip->dbx== -1) { printf("randspk() DBXA: %d %g %g",ip->rvi,next,WEX); chk(2.); }
  }
  ENDVERBATIM
  
}


PROCEDURE vers () {
  printf("$Id
}


VERBATIM
void val(double xx, double ta) {
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  if (VGBdel>0) {
    vii[4]=esinr(xx-tGB);
  } else {
    vii[4]=VGB*EXP(-(xx - ta)/tauGB);
  }  
  vii[5]=AHP*EXP(-(xx - ta)/tauahp);
  vii[6]=vii[1]+vii[2]+vii[3]+vii[4]+vii[5];
}
ENDVERBATIM


VERBATIM
void valps(double xx, double ta) {
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  
  
  vii[6]=vii[1]+vii[2]+vii[3]; 
}
ENDVERBATIM


PROCEDURE record () {
  VERBATIM {
  int k; double ti;
  vp = SOP;
  if (tg>=t) return 0;
  if (vp->p >= vp->size) {
    if (errflag) return 0;
    printf("**** WARNING out of recording room for INTF type%d id%d at %g****\n",IDP->type,IDP->id,t);
    printf("**************** WARNING: No further WARNINGS ****************\n");
    errflag=1;
    return 0;
  }
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
  if (vp->p >= vp->size || vp->vvo[6]==0) {
    return 0;
  }
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
  {int max, i,err;
    double last;
    if (! ifarg(1)) {printf("Return initvspks(ivspks,vspks,wvspks)\n"); return 0.;}
    ip=IDP;  err=0;
    i = vector_arg_px(1, &isp); 
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("initvspks ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("initvspks ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("initvspks ERR: 3rd vec is of different size\n");}
    vspn=max;
    ip->vinflg=1;
    for (i=0; i<max && (int)isp[i] != ip->id ; i++); 
    if (i==max) { 
      printf("initvspks WARN: %d not found in ivspks\n",ip->id); 
      ip->vinflg=0; ip->rve=-1;
      return(0.); 
    }
    ip->rvb=ip->rvi=i;
    last=vsp[i++];
    for (; i<max && (int)isp[i] == ip->id ; i++) { 
      if (vsp[i]<=last) { err=1; 
        printf("initvspks ERR: nonmonotonic for cell#%d: %g %g\n",ip->id,last,vsp[i]); }
      last=vsp[i];
    }
    ip->rve=i-1;
    if (subsvint>0) { 
      vsp[ip->rve] = vsp[ip->rvb]+subsvint;
      wsp[ip->rve] = wsp[ip->rvb];
    }
    if (err) { ip->rve=0; hoc_execerror("",0); }
  }
  ENDVERBATIM
}




PROCEDURE trvsp ()
{
  VERBATIM 
  int i, flag; 
  double ind, t0_local;
  ip=IDP;
  flag=(int) *getarg(1);
  if (subsvint==0.) {printf("trvsp"); return(0.);}
  ind=isp[0]; t0_local=vsp[0];
  if (flag==1) {
    for (i=0; i<vspn; i++) {
      if (isp[i]!=ind) {
        vsp[i-1]=1.e9;
        ind=isp[i];
      }
    }
    vsp[vspn-1]=1.e9;
  } else if (flag==2) {
    for (i=0; i<vspn; i++) {
      if (isp[i]!=ind) {
        vsp[i-1] = t0_local + subsvint;
        ind = isp[i];
        t0_local = vsp[i];
      }
    }
    vsp[vspn-1] = t0_local + subsvint;
  } else {printf("trvsp flag %d not recognized\n",flag); hxe();}
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
  printf("initinvl() NOT BEING USED\n")
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
  if (IDP->wrec) { wrecord(1e9); }
  if (IDP->record) {
    record(); 
    for (k=0;k<NSV;k++) if (vp->vvo[k]!=0) { 
      vector_resize(vp->vv[k], vp->p);
    }
  }}
  ENDVERBATIM
}



PROCEDURE chk (f) {
  VERBATIM 
  {int i,lfg;
  lfg=(int)_lf;
  ip=IDP;
  printf("ID:%d; typ: %d; rec:%d wrec:%d inp:%d jit:%d invl:%d\n",ip->id,ip->type,ip->record,ip->wrec,ip->input,ip->jitter,ip->invl0);
  if (lfg==1) {
    if (SOP!=nil) {
      vp = SOP;
      printf("p %d size %d tg %g\n",vp->p,vp->size,tg);
      for (i=0;i<NSV;i++) printf("%d %x %x;",i,vp->vv[i],vp->vvo[i]);
    } else printf("Recording pointers not initialized");
  }
  if (lfg==2) { 
    printf("Global vectors for input and jitter: \n");
    if (vsp!=nil) printf("VSP: %x (%d/%d-%d)\n",vsp,ip->rvi,ip->rvb,ip->rve); else printf("no VSP\n");
    if (jsp!=nil) printf("JSP: %x (%d/%d)\n",jsp,jtpt,jitmax); else printf("no JSP\n");
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
    printf("wwwid %g wwht %d nsw %g\n WW vecs: ",wwwid,(int)wwht,nsw);
    for (i=0;i<NSW;i++) printf("%d %x %x;",i,ww[i],wwo[i]);
  }}
  ENDVERBATIM
}


FUNCTION pid () {
  VERBATIM 
  printf("INTF%d(%d/%d@%g) ",IDP->id,IDP->type,IDP->col,t);
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

FUNCTION col () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->col = (unsigned char) *getarg(1);
  _lcol = (double)ip->col;
  ENDVERBATIM
}

FUNCTION dbx () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->dbx = (unsigned char) *getarg(1);
  _ldbx = (double)ip->dbx;
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



PROCEDURE initwrecOLD () {
  VERBATIM 
  {int k;
  if (! ifarg(4)) { 
    wwsz=WSZ;
    wf1 = hoc_obj_file_arg(1);
    wf2 = hoc_obj_file_arg(2);
  } else if (! ifarg(8)) { 
    if (WSW!=4) { 
      printf("INTF initwrec ERR w-vecs compiled for 4 args\n");
      hoc_execerror("",0); }
    IDP->wrec=1;
    for (k=0;k<WSW;k++) {
      ww[k]=vector_arg(k+1);
      wwaz=vector_arg_px(k+1, &(wwo[k]));
    }
    if (wwsz==0) wwsz=(unsigned int)vector_buffer_size(ww[0]); 
    for (k=0;k<WSW;k++) if (wwsz!=(unsigned int)vector_buffer_size(ww[k])) {
      printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,vector_buffer_size(ww[k]));
    }
  } else { 
    if (WSW!=8) { 
      printf("INTF initwrec ERR w-vecs compiled for 8 args\n");
      hoc_execerror("",0); }
    IDP->wrec=1;
    for (k=0;k<WSW;k++) {
      ww[k]=vector_arg(k+1);
      wwaz=vector_arg_px(k+1, &(wwo[k]));
    }
    if (wwsz==0) wwsz=(unsigned int)vector_buffer_size(ww[0]); 
    for (k=0;k<WSW;k++) if (wwsz!=(unsigned int)vector_buffer_size(ww[k])) {
      printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,vector_buffer_size(ww[k]));
    }
  }}
  ENDVERBATIM
}


PROCEDURE initwrec () {
  VERBATIM 
  {int i, k, num, cap;  Object* ob;
    ob =   *hoc_objgetarg(1); 
    num = ivoc_list_count(ob);
    if (num>NSW) { printf("INTF initwrec() WARN: can only store %d ww vecs\n",NSW); hxe();}
    nsw=(double)num;
    for (k=0;k<num;k++) {
      cap = list_vector_px2(ob, k, &wwo[k], &ww[k]);
      if (k==0) wwsz=cap; else if (wwsz!=cap) {
        printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,cap); hxe(); }
    }
  }
  ENDVERBATIM
}


PROCEDURE wrecordOLD (t,s0,w0) {
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



PROCEDURE popspk (x) {
  TABLE Psk DEPEND wwwid,wwht FROM -40 TO 40 WITH 81
  Psk = -wwht*exp(-2.*x*x/wwwid/wwwid)
}

PROCEDURE pskshowtable () {
  VERBATIM 
  int j;
  printf("_tmin_popspk:%g -_tmin_popspk:%g\n",_tmin_popspk,-_tmin_popspk);
  for (j=0;j<=-2*(int)_tmin_popspk+1;j++) printf("%g ",_t_Psk[j]);
  printf("\n");
  ENDVERBATIM 
}


PROCEDURE wrecord (te) {
  VERBATIM 
  {int j,k,max,wrp; double ti,scale;
  wrp=(int)IDP->wrec-1; 
  scale=(double)IDP->wscale;
  if (_lte<1.e9) { 
    max=-(int)_tmin_popspk; 
    k=(int)floor((_lte-rebeg)/vdt+0.5);
    for (j= -max;j<=max && k+j>0 && k+j<wwsz;j++) {
      wwo[wrp][k+j] += scale*_t_Psk[j+max]; 
    }
  } else if (twg>=t) {
    return 0;
  } else {
    for (ti=twg,k=(int)floor((twg-rebeg)/vdt+0.5);ti<=t && k<wwsz;ti+=vdt,k++) { 
      valps(ti,twg);  
      wwo[wrp][k]+=scale*vii[6];
    }
    twg=ti;
  }
  }
  ENDVERBATIM
}

FUNCTION wrec () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->wrec = (unsigned char) *getarg(1);
  if (ifarg(2)) ip->wscale = (float) *getarg(2); else ip->wscale=1.;
  _lwrec=(double)ip->wrec;
  ENDVERBATIM
}

FUNCTION wwszset () {
  VERBATIM
  if (ifarg(1)) wwsz = (unsigned int) *getarg(1);
  _lwwszset=(double)wwsz;
  ENDVERBATIM
}


FUNCTION wwfree () {
  VERBATIM
  int k;
  IDP->wrec=0;
  wwsz=0; wwpt=0; nsw=0.;
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
  popspk(0) 
  VERBATIM
  int j,k;
  if (nsw>0. && wwo[0]!=0) { 
    printf("Initializing ww to record for %g (%g)\n",vdt*wwsz,vdt);
    wwpt=0;
    for (k=0;k<(int)nsw;k++) {
      vector_resize(ww[k], wwsz);
      for (j=0;j<wwsz;j++) wwo[k][j]=0.;
    }
  } else printf("global_init WARNING: wrec not initialized\n");
  ENDVERBATIM
}

PROCEDURE global_fini () {
  VERBATIM
  int k;
  for (k=0;k<(int)nsw;k++) vector_resize(ww[k], (int)floor(t/vdt+0.5));
  ENDVERBATIM
}

PROCEDURE global_finiOLD () {
  VERBATIM
  {int k;
  if (IDP->wrec) {
    if (wwo[0]!=0) { 
      for (k=0;k<WSW;k++) vector_resize(ww[k], wwpt);
    } else {
      fprintf(wf1,"
      fwrite(&wwt,sizeof(float),wwpt,wf2);  
      fwrite(&wwi,sizeof(int),wwpt,wf2);  
      fwrite(&wws,sizeof(char),wwpt,wf2);  
      fwrite(&www,sizeof(float),wwpt,wf2);  
      printf("Closing file with wwpt=%d at location %d\n",wwpt,ftell(wf2));
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
  if (ifarg(1)) ip->jitter = (unsigned char) *getarg(1);
  _ljitset=(double)ip->jitter;
  ENDVERBATIM
}



FUNCTION flag () {
  VERBATIM
  char *sf; int ii;
  ip = IDP;
  sf = gargstr(1);
  for (ii=0;ii<iflnum && strncmp(sf, &iflags[ii*4], 3)!=0;ii++) ;
  if (ii==10) {printf("INTF ERR: %s not found as a flag (%s)\n",sf,iflags); hxe();}
  if (ifarg(2)) (&ip->type)[ii] = (unsigned char) *getarg(2);  
  _lflag=(double)(unsigned char)(&ip->type)[ii];
  ENDVERBATIM
}


FUNCTION invlset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->invl0 = (unsigned char) *getarg(1);
  _linvlset=(double)ip->invl0;
  ENDVERBATIM
}


FUNCTION vinset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->vinflg = (unsigned char) *getarg(1);
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