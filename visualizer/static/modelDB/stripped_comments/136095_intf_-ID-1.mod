NEURON {
  ARTIFICIAL_CELL INTF
  RANGE VAM, VNM, VGA, VGB, AHP           
  RANGE Vm                                
  
  RANGE tauAM, tauNM, tauGA, tauGB     
  RANGE tauahp, ahpwt                  
  RANGE tauRR , RRWght                 
  RANGE RMP,VTH,Vblock,VTHC,VTHR       
  RANGE nbur,tbur,refrac               
  RANGE invl,oinvl,WINV,invlt           
  RANGE VGBdel,tGB,VGBa,Vbrefrac,rebound,rebob,offsetGB   
  RANGE STDAM, STDNM, STDGA             
  GLOBAL EAM, ENM, EGA, EGB, mg         
  GLOBAL tauGBGP,wGBGP,GPkd,Gn          
  GLOBAL spkht, wwwid,wwht              
  GLOBAL stopoq                         
  
  POINTER sop                          
  RANGE  spck,xloc,yloc,zloc
  RANGE  t0,tg,twg,tGB,refractory,trrs 
  RANGE  cbur,WEX                      
  GLOBAL vdt,nxt,RES,ESIN,Bb,Psk      
  GLOBAL prnum, nsw, rebeg             
  GLOBAL subsvint, jrsvn, jrsvd, jrtime, jrtm 
  GLOBAL DEAD_DIV, seedstep            
  GLOBAL seaddvioff                    
  GLOBAL WVAR,DELMIN
  GLOBAL savclock,slowset,FLAG  
  GLOBAL tmax,installed,verbose        
  GLOBAL pathbeg,pathend,PATHMEASURE,pathidtarg,pathtytarg,seadsetting,pathlen
}

PARAMETER {
  tauAM = 10 (ms)
  tauNM = 300 (ms)
  tauGA = 10 (ms)
  tauGB = 300 (ms)
  tauGBGP = 50 (ms) 
  invl =  100 (ms)
  WINV =  0
  wGBGP = 1 (ms) 
  GPkd  = 100    
  ahpwt = 0
  tauahp= 10 (ms)
  tauRR = 6 (ms)
  refrac = 5 (ms)
  Vbrefrac = 20 (ms)
  RRWght = 0.75
  wwwid = 10
  wwht = 10
  VTH = -45      
  VTHC = -45
  VTHR = -45
  Vblock = -20   
  vdt = 0.1      
  mg = 1         
  sop=0
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
  rebeg=0
  subsvint=0
  jrsvn=1e4 jrsvd=1e4 jrtime=-1 jrtm=-1
  seedstep=44340
  seaddvioff=9102098713763e-134
  DEAD_DIV=1
  WVAR=0.2
  stopoq=0
  PATHMEASURE=0
  verbose=1
  seadsetting=0
  pathidtarg=-1
  DELMIN=1e-5 
  STDAM=0
  STDNM=0
  STDGA=0
}

ASSIGNED {
  Vm VAM VNM VGA VGB AHP
  VGBa t0 tGB tg twg refractory nxt xloc yloc zloc trrs
  WEX RES ESIN Gn Bb Psk cbur invlt oinvl rebob tmax spck savclock slowset FLAG
  installed
  pathbeg pathend pathtytarg pathlen
}



CONSTRUCTOR {
  VERBATIM 
  { int lid,lty,lin,lco,i; unsigned int sz;
    if (ifarg(1)) { lid=(int) *getarg(1); } else { lid= UINT_MAX; }
    if (ifarg(2)) { lty=(int) *getarg(2); } else { lty= -1; }
    if (ifarg(3)) { lin=(int) *getarg(3); } else { lin= -1; }
    if (ifarg(4)) { lco=(int) *getarg(4); } else { lco= -1; }
    _p_sop = (double*)ecalloc(1, sizeof(id0)); 
    ip = IDP;
    ip->id=lid; ip->type=lty; ip->inhib=lin; ip->col=lco; 
    ip->pg=0x0; ip->dvi=0x0; ip->sprob=0x0;
    ip->dead = ip->invl0 = ip->record = ip->jttr = ip->input = 0; 
    ip->dvt = ip->vbr = ip->wrec = ip->jcn = 0;
    for (i=0;i<WRNUM;i++) {ip->wreci[i]=-1; ip->wscale[i]=-1.0;}
    ip->rve=-1;
    pathbeg=-1;
    slowset=jrmax=0; 
    process=(int)getpid();
    CNAME[SU]="SU"; CNAME[DP]="DP"; CNAME[IN]="IN";
    if (installed==2.0) { 
      sz=ivoc_list_count(ce);
      printf("\t**** WARNING new INTF created: may want to rerun jitcondiv ****\n");
    } else installed=1.0; 
    cbsv=0x0;
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
  free(IDP);
  }
  ENDVERBATIM
}


INITIAL { LOCAL id
  reset() 
  t0 = 0
  tg = 0
  twg = 0
  tGB = 0
  trrs = 0
  tmax=0
  pathend=-1
  pathlen=0
  spck=0
  VERBATIM
  { int i,ix;
  ip=IDP;
  _lid=(double)ip->id;
  ip->spkcnt=0;
  ip->blkcnt=0;
  ip->errflag=0;
  for (i=0;i<CTYN;i++){ix=cty[i]; blockcnt[ix]=spikes[ix]=AMo[ix]=NMo[ix]=GAo[ix]=GBo[ix]=0;}
  }
  ENDVERBATIM
  jrsvn=jrsvd jrtime=jrtm
  
  if (vinflag()) { randspk() net_send(nxt,2)}
  if (recflag()) { recini() } 
  if (pathbeg==id) { 
    stoprun=0 
    net_send(0,2) 
  } 
  rebeg=0 
}

PROCEDURE reset () {
  Vm = RMP
  VAM = 0
  VNM = 0
  VGA = 0
  VGB = 0
  VGBa = 0
  offsetGB=0
  AHP=0
  rebob=-1e9
  invlt = -1
  t0 = t
  tGB = t
  tg = t
  twg = t
  trrs = t
  cbur = 0 
  spck = 0 
  refractory = 0 
  VTHC=VTH 
  VTHR=VTH 
}

VERBATIM
unsigned int GetDVIDSeedVal(unsigned int id) {
  double x[2];
  if (seadsetting==1) { 
    sead=((unsigned int)ip->id+seaddvioff)*1e6;
  } else { 
    if (seadsetting==2) printf("Warning: GetDVIDSeedVal called with wt rand turned off\n");
    x[0]=(double)id; x[1]=seaddvioff;
    sead=hashseed2(2, x);
  }
  return sead;
}
ENDVERBATIM


FUNCTION DVIDSeed(){
  VERBATIM
  return (double)GetDVIDSeedVal(IDP->id);
  ENDVERBATIM
}


NET_RECEIVE (wAM,wNM,wGA,wGB,wflg) { LOCAL tmp,jcn,id
  INITIAL { wNM=wNM wGA=wGA wGB=wGB wflg=0}
  
VERBATIM
  int prty,poty,prin,ii,sy,nsyn; double STDf; 

ENDVERBATIM
  tmax=t
  VERBATIM
  if (stopoq && !qsz()) stoprun=1;
  ip=IDP; pg=ip->pg;
  if (ip->dead) return; 
  _ljcn=ip->jcn; _lid=ip->id;
  tpnt = _pnt; 
  if (PATHMEASURE) { 
    if (_lflag==2 || _lflag<0) { 
      double idty; int i;
      if (_lflag==2) ip->flag=-1; 
      idty=(double)(FOFFSET+ip->id)+1e-2*(double)ip->type+1e-3*(double)ip->inhib+1e-4;
      for (i=0;i<ip->dvt && !stoprun;i++) if (ip->sprob[i]) {
        (*pnt_receive[get_type(ip->dvi[i]->_prop)])(ip->dvi[i], wts, idty);
        
#ifdef NRN_MECHANISM_DATA_IS_SOA
        neuron::legacy::set_globals_from_prop(_pnt->_prop, _ml_real, _ml, _iml);
#else
        _p = _pnt->_prop->param;
#endif
        _ppvar = get_dparam(_pnt->_prop);
        ip = IDP;
      }
      return;  
    } else if (_lflag!=2 && (pathtytarg==(double)ip->type || pathidtarg==(double)ip->id)) {
      if (pathend==(double)ip->id) return; 
      ip->flag=(unsigned char)floor(t)+1; 
      pathend=(double)ip->id; 
      pathlen=tmax+1; 
      stoprun=1.; 
      return;
      
    } else if (ip->flag   || ip->dvt==0 || stoprun) {
      return; 
    } else if (ip->inhib) {
      if (!ip->flag) ip->flag=(unsigned char)floor(t)+1;
    } else { 
      ip->flag=(unsigned char)floor(t)+1;
   #if defined(t)
      net_send((void**)0x0, wts,tpnt,t+1.,-1.); 
  #else
      net_send((void**)0x0, wts,tpnt,1.,-1.); 
  #endif
      return;
    }
  }
  if (_lflag==OK) { FLAG=OK; flag(); return; } 
  if (_lflag<0) { callback(_lflag); return; }
  eventtot+=1;
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


  if (flag>=FOFFSET) { 
    VERBATIM {
      
      poty=(int)ip->type;
      prty=(int)(1e2*(_lflag-floor(_lflag)));
      prin=(int)(1e3*(_lflag-floor(_lflag)-prty*1e-2)); 
      sy=prin?GA:AM;
      STDf=_args[0]; 
      for (ii=0;ii<=4;ii++) _args[ii]=0.; 
      for (ii=sy,nsyn=0;ii<sy+2;ii++) nsyn+=((_args[ii]=WMAT(prty,poty,ii)*WD0(prty,poty,ii))>0.);
      if (nsyn==0) return;
      if (seadsetting!=2) { 
        if (seadsetting==1) {
          sead=(unsigned int)(floor(_lflag)*ip->id*seedstep); 
        } else { 
          hsh[0]=floor(_lflag); hsh[1]=(double)ip->id; hsh[2]=seedstep;
          sead=hashseed2(3, hsh); 
        }
        mcell_ran4(&sead, &_args[sy], 2, 1.);
        for (ii=sy;ii<sy+2;ii++) { 
          _args[ii]=2*WVAR*(_args[ii]+0.5/WVAR-0.5)*WMAT(prty,poty,ii)*WD0(prty,poty,ii);
        }
      }
    }
    ENDVERBATIM
VERBATIM
    if (ip->dbx>2) 
ENDVERBATIM
{ 
      pid() 
      printf("DF
      printf(" (%g %g %g %g %g)",wAM,wNM,wGA,wGB,wflg)
      printf("\n")
    }

  } else if (flag==4) { 
    cbur=cbur-1  
    if (cbur>0) { 
      net_send(tbur,4) 
    } else { 
      refractory = 1      
      net_send(refrac, 3) 
    }
    tmp=t
VERBATIM
    if (ip->jttr) 
ENDVERBATIM
{ tmp= t+jttr()/10 } 
    if (jcn) { jitcon(tmp) } else { net_event(tmp) }
VERBATIM
    spikes[ip->type]++; 

ENDVERBATIM
    spck=spck+1
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
VERBATIM
    return; 

ENDVERBATIM
    
    
    
    
  } else if (flag==0 && wGB==0 && wflg==1) {
VERBATIM
    ip->input=1; 

ENDVERBATIM
    wflg=2 
    randspk() 
    net_send(nxt,2)
VERBATIM
    return; 

ENDVERBATIM
  } else if (flag==0 && wGB==0 && wflg==2) { 
VERBATIM
    ip->input=0; 

ENDVERBATIM
    wflg=1  
VERBATIM
    return; 

ENDVERBATIM
  }
  
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
  if(refractory==0){
    if(VTHC>VTH) { VTHC = VTH + (VTHR-VTH)*EXP(-(t-trrs)/tauRR) }
  }
  if (VGBdel>0) {
    VGB = esinr(t-tGB) 
  } else {
    if (VGB< -hoc_epsilon){ 
      VGB = VGB*EXP(-(t - t0)/tauGB) } else { VGB=0 }
  }      
  if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } 
  t0 = t 
  Vm = VAM+VNM+VGA+VGB+AHP 
  if (Vm> -RMP) {Vm= -RMP}
  if (Vm<  RMP) {Vm= RMP} 

  if (flag==0 || flag>=FOFFSET) { 
    
    if (wAM>0) {
      if (rebob!=1e9 && rebob!=-1e9) {
VERBATIM
        cbur=floor(rebound*rebob/EGB); 

ENDVERBATIM
        net_send(tbur,4) 
        rebob=1e9
      }
      if (STDAM==0) { VAM = VAM + wAM*(1-Vm/EAM)
      } else        { VAM = VAM + (1-STDAM*STDf)*wAM*(1-Vm/EAM) }
      if (VAM>EAM) { 
VERBATIM
        AMo[ip->type]++; 

ENDVERBATIM
      } else if (VAM<0) { VAM=0 }
    }
    
    if (wNM>0 && VNM<ENM) { 
      rates(RMP+Vm)
      if (STDNM==0) { VNM = VNM + wNM*Bb*(1-Vm/ENM) 
      } else        { VNM = VNM + (1-STDNM*STDf)*wNM*Bb*(1-Vm/ENM) }
      if (VNM>ENM) { 
VERBATIM
        NMo[ip->type]++; 

ENDVERBATIM
      } else if (VNM<0) { VNM=0 }
    }
    
    if (wGA>0 && VGA>EGA) { 
      if (STDGA==0) {  VGA = VGA - wGA*(1-Vm/EGA) 
      } else {         VGA = VGA - (1-STDGA*STDf)*wGA*(1-Vm/EGA) }
      if (VGA<EGA) { 
VERBATIM
        GAo[ip->type]++; 

ENDVERBATIM
VERBATIM
  if (ip->dbx>2) 
ENDVERBATIM
{ 
    pid() printf("DB0A
    if (flag==0) { printf(" (%g %g %g %g %g %g)",wGA,EGA,VGA,Vm,AHP,STDf) }  
VERBATIM
    printf("\nAA:%d:%d\n",GAo[ip->type],ip->type); 

ENDVERBATIM
    printf("\n")
  }
      } else if (VGA>0) { VGA=0 } 
    }
    if (wGB>1e-6) {
      if (VGBdel>0) { net_send(VGBdel,5)  
      } else { 
        
        wflg=wflg*EXP(-(t-tGB)/tauGBGP)+wGBGP 
        coop(wflg)               
        VGB = VGB - wGB*(1-Vm/EGB)*Gn
        if (VGB<EGB) { 
VERBATIM
          GBo[ip->type]++; 

ENDVERBATIM
        }
        if (VGB<rebob && rebob!=1e9 && rebob!=-1e9) { rebob=VGB }
      }
    }

VERBATIM
    if (ip->invl0) 
ENDVERBATIM
{ 
      Vm = RMP+VAM+VNM+VGA+VGB+AHP
      if (Vm>0)   {Vm= 0 }
      if (Vm<-90) {Vm=-90}
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
  } else if (flag==1) { 
    
    if (WINV<0) { 
      if (jcn) { jitcon(t) } else { net_event(t) } 
VERBATIM
      spikes[ip->type]++; 

ENDVERBATIM
      spck=spck+1
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
    if (flag==2) { 
{pid() printf("DBBa
      if (WEX<0) { 
        if (jcn) { jitcon(t) } else { net_event(t) } 
VERBATIM
        spikes[ip->type]++; 

ENDVERBATIM
        spck=spck+1
VERBATIM
        if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBB
        if (WEX<-1 && WEX!=-1e9) { cbur=-WEX  net_send(tbur,4) }
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
      if (WEX!=-1e9) { 
        randspk() 
        if (nxt>0) { net_send(nxt,2) }
      }
    }
  } else if (flag==3) { 
    refractory = 0 
    trrs = t 
VERBATIM
    return; 

ENDVERBATIM
  }

  Vm = VAM+VNM+VGA+VGB+RMP+AHP 
  if (Vm>0)   {Vm= 0 }
  if (Vm<-90) {Vm=-90}
  if (refractory==0 && Vm>VTHC) {
VERBATIM
    if (!ip->vbr && Vm>Vblock) {

ENDVERBATIM
VERBATIM
      ip->blkcnt++; blockcnt[ip->type]++; return; }

ENDVERBATIM
    AHP = AHP - ahpwt
    tmp=t
    
VERBATIM
    if (ip->jttr) 
ENDVERBATIM
{ tmp= t+jttr() }  
    if (jcn) { jitcon(tmp) } else { net_event(tmp) } 
VERBATIM
    spikes[ip->type]++; 

ENDVERBATIM
    spck=spck+1
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
    VTHC=VTH+RRWght*(Vblock-VTH)
    VTHR=VTHC 
    refractory = 1 
    if (nbur>1) { 
      cbur=nbur-1 net_send(tbur,4) 
VERBATIM
      return; 

ENDVERBATIM
    } else if (rebob==1e9) { rebob=0 }
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
    net_send(refrac, 3) 
  }
}








PROCEDURE jitcon (tm) {
  VERBATIM {
  double mindel, randel, idty, *x; int prty, poty, i, j, k, dv; 
  Point_process *pnt; void* voi;
  
  
  ip=IDP;
  if (!pg) {printf("No network defined -- must run jitcondiv() and pgset\n"); hxe();}
  ip->spkcnt++; 
  if (jrj<jrmax) { 
    jrid[jrj]=(double)ip->id; jrtv[jrj]=_ltm;
    jrj++;
  } else if (wf2 && jrmax) spkoutf2(); 
  jri++;  
  if (jrtm>0) {
    if (t>jrtime) {
      jrtime+=jrtm;
      spkstats2(1.);
    }
  } else if (jrsvd>0 && jri>jrsvn) { 
    jrsvn+=jrsvd; printf("t=%.02f %ld ",t,jri);
    spkstats2(1.);
  }
  prty=(int)ip->type;
  if (ip->jcn==1) if (ip->dvt>0) {  
      #if defined(t)
    if (ip->jcn==1) if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,t+ip->del[0],-1.);
      #else
    if (ip->jcn==1) if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,ip->del[0],-1.);
      #endif
  }
  }   
  ENDVERBATIM  
}


PROCEDURE spkstats () {
VERBATIM {
  if (ifarg(1)) tf=hoc_obj_file_arg(1); else tf=0x0;
}
ENDVERBATIM
}


PROCEDURE spkoutf () {
VERBATIM {
  if (ifarg(2)) {
    wf1=hoc_obj_file_arg(1); 
    wf2=hoc_obj_file_arg(2);
  } else if (wf1 != 0x0) {
    spkoutf2();
    wf1=(FILE*)0x0; wf2=(FILE*)0x0;
  }
}
ENDVERBATIM
}

VERBATIM
static void spkoutf2 () {
    fprintf(wf1,"
    fwrite(jrtv,sizeof(double),jrj,wf2); 
    fwrite(jrid,sizeof(double),jrj,wf2); 
    fflush(wf1); fflush(wf2);
    jrj=0;
}
ENDVERBATIM

PROCEDURE callhoc () {
  VERBATIM
  if (ifarg(1)) {
    cbsv=hoc_lookup(gargstr(1));
  } else {
    cbsv=0x0;
  }
  ENDVERBATIM
}


PROCEDURE spkstats2 (flag) {
VERBATIM {
  int i, ix, flag; double clk;
  flag=(int)(_lflag+1e-6);
  clk=clock()-savclock; savclock=clock();
  if (cbsv) hoc_call_func(cbsv,0);
  if (tf) fprintf(tf,"t=%.02f;%ld(%g) ",t,jri,clk/1e6); else {
    printf("t=%.02f;%ld(%g) ",t,jri,clk/1e6); }
  for (i=0;i<CTYN;i++) {
    ix=cty[i];
    spktot+=spikes[ix];
    if (tf) {
      fprintf(tf,"%s:%d/%d:%d;%d;%d;%d ",CNAME[i],spikes[ix],\
              blockcnt[ix],AMo[ix],NMo[ix],GAo[ix],GBo[ix]);
    } else {
      printf("%s:%d/%d:%d;%d;%d;%d ",CNAME[i],spikes[ix],blockcnt[ix],\
             AMo[ix],NMo[ix],GAo[ix],GBo[ix]);
    }
    spck=0;
    blockcnt[ix]=spikes[ix]=0;
    AMo[ix]=NMo[ix]=GAo[ix]=GBo[ix]=0;
  }
  if (tf && flag==2) {  fprintf(tf,"\nt=%g tot_spks: %ld; tot_events: %ld\n",t,spktot,eventtot); 
  } else if (flag==2) {  printf("\ntotal spikes: %ld; total events: %ld\n",spktot,eventtot); 
  } else if (tf) fprintf(tf,"\n"); else printf("\n");
}
ENDVERBATIM
}

PROCEDURE oobpr () {
VERBATIM {
  int i,ix;
  for (i=0;i<CTYN;i++){ 
    ix=cty[i];
    printf("%d:%d/%d:%d;%d;%d;%d ",ix,spikes[ix],blockcnt[ix],AMo[ix],NMo[ix],GAo[ix],GBo[ix]);
  }
  printf("\n");
}
ENDVERBATIM
}

PROCEDURE callback (fl) {
  VERBATIM {
  int i; double idty, del0, ddel; id0 *jp; Point_process *upnt; 
  i=(unsigned int)((-_lfl)-1); 
  jp=IDP; upnt=tpnt; del0=jp->del[i]; ddel=0.;
  idty=(double)(FOFFSET+jp->id)+1e-2*(double)jp->type+1e-3*(double)jp->inhib+1e-4;
  while (ddel<=DELMIN) { 
    if (Vblock<VTHC) { 
      wts[0]=0; 
    } else { 
      wts[0]=(VTHC-VTH)/(Vblock-VTH); 
    }
    if (jp->sprob[i]) (*pnt_receive[get_type(jp->dvi[i]->_prop)])(jp->dvi[i], wts, idty);
    
#ifdef NRN_MECHANISM_DATA_IS_SOA
    neuron::legacy::set_globals_from_prop(upnt->_prop, _ml_real, _ml, _iml);
#else
    _p = upnt->_prop->param;
#endif
    _ppvar = get_dparam(upnt->_prop);
    i++;
    if (i>=jp->dvt) return 0; 
    ddel=jp->del[i]-del0;   
  }
  
  while (i<jp->dvt && (!jp->sprob[i] || id0ptr(jp->dvi[i]->_prop)->dead)) i++;
  if (i<jp->dvt) {
    ddel= jp->del[i] - del0;;
  #if defined(t)
    net_send((void**)0x0, wts,upnt,t+ddel,(double) -(i+1)); 
  #else
    net_send((void**)0x0, wts,upnt,ddel,(double) -(i+1)); 
  #endif
  }
  } 
  ENDVERBATIM
}



PROCEDURE mkdvi () {
VERBATIM {
  int i,j,k,prty,poty,dv,dvt,dvii; double *x, *db, *dbs; 
  Object *lb;  Point_process *pnnt, **da, **das;
  ip=IDP; pg=ip->pg; 
  if (ip->dead) return 0;
  prty=ip->type;
  sead=GetDVIDSeedVal(ip->id);
  for (i=0,k=0,dvt=0;i<CTYN;i++) { 
    poty=cty[i];
    dvt+=DVG(prty,poty);
  }
  da =(Point_process **)malloc(dvt*sizeof(Point_process *));
  das=(Point_process **)malloc(dvt*sizeof(Point_process *)); 
  db =(double *)malloc(dvt*sizeof(double)); 
  dbs=(double *)malloc(dvt*sizeof(double)); 
  for (i=0,k=0,dvii=0;i<CTYN;i++) { 
    poty=cty[i];
    dv=DVG(prty,poty);
    if (dv>0) {
      sead+=dv;
      if (dv>dscrsz) {
        printf("B:Divergence exceeds dscrsz: %d>%d for %d->%d\n",dv,dscrsz,prty,poty); hxe(); }
      mcell_ran4(&sead, dscr ,  dv, pg->ixe[poty]-pg->ix[poty]+1);
      for (j=0;j<dv;j++) {
        if (!(lb=ivoc_list_item(ce,(unsigned int)floor(dscr[j]+pg->ix[poty])))) {
          printf("INTF:callback %g exceeds %d for list ce\n",floor(dscr[j]+pg->ix[poty]),cesz); 
          hxe(); }
        pnnt=(Point_process *)lb->u.this_pointer;
        da[j+dvii]=pnnt;
      }
      mcell_ran4(&sead, dscr , dv, 2*DELD(prty,poty));
      for (j=0;j<dv;j++) {
        db[j+dvii]=dscr[j]+DELM(prty,poty)-DELD(prty,poty); 
        if (db[j+dvii]<0) db[j+dvii]=-db[j+dvii];
      }
      dvii+=dv;
    }
  }
  gsort2(db,da,dvt,dbs,das);
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt;
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); 
  for (i=0;i<dvt;i++) ip->sprob[i]=1; 
  free(da); free(db); 
  }
ENDVERBATIM
}


PROCEDURE patha2b () {
  VERBATIM
  int i; double idty, *x; static Point_process *_pnt; static id0 *ip0;
  pathbeg=*getarg(1); pathidtarg=*getarg(2);
  pathtytarg=-1;  PATHMEASURE=1; pathlen=stopoq=0;
  for (i=0;i<cesz;i++) { lop(ce,i); 
    if ((i==pathbeg || i==pathidtarg) && qp->inhib) {
      pid(); printf("Checking to or from inhib cell\n" ); hxe(); }
    qp->flag=qp->vinflg=0; 
  }
  hoc_call_func(hoc_lookup("finitialize"), 0);
  cvode_fadvance(1000.0); 
  ENDVERBATIM
}



FUNCTION pathgrps () {
  VERBATIM
  int i,j,k,na,nb,flag; double idty,*a,*b,*x,sum; static Point_process *_pnt; static id0 *ip0;
  Symbol* s; char **pfl;
  x=0x0;
  s=hoc_lookup("finitialize");
  if (ifarg(2)) {
    na=vector_arg_px(1,&a);
    nb=vector_arg_px(2,&b);
    if (ifarg(3)) x=vector_newsize(vector_arg(3),na*nb);
  } else {
    na=nb=cesz;  
    if (ifarg(1)) x=vector_newsize(vector_arg(1),na*nb);
  }
  
  pfl = (char **)malloc(cesz * (unsigned)sizeof(char *));
  for (i=0;i<cesz;i++) { lop(ce,i); scr[i]=qp->inhib; pfl[i]=&qp->flag; }
  pathtytarg=-1;  PATHMEASURE=1; pathlen=stopoq=0;
  for (k=0,sum=0;k<na;k++) {
    pathbeg=a[k]; 
    if (scr[(int)pathbeg]) { 
      if (x) for (j=0;j<nb;j++) x[k*nb+j]=0.;
      continue;
    }
    for (j=0;j<nb;j++) { 
      pathidtarg=b[j]; 
      if (scr[(int)pathidtarg]) { if (x) x[k*nb+j]=0.; 
        continue;
      }
      
      for (i=0;i<cesz;i++) *pfl[i]=0;
      hoc_call_func(s, 0);
      cvode_fadvance(1000.0); 
      sum+=pathlen;
      if (x) x[k*nb+j]=pathlen;
    }
  }
  PATHMEASURE=0;
  free(pfl);
  _lpathgrps=sum/na/nb;
  ENDVERBATIM
}







FUNCTION getdvi () {
  VERBATIM 
  {
    int i,j,k,iarg,av1,a2,a3,a4,dvt,getactive=0,idx=0,*pact,prty,poty,sy,ii; 
    double *dbs, *x,*x1,*x2,*x3,*x4,*x5,idty,y[2],flag;
    void* voi, *voi2,*voi3; Point_process **das;
    ip=IDP; pg=ip->pg; 
    getactive=a2=a3=a4=0;
    if (ip->dead) return 0;
    dvt=ip->dvt;
    dbs=ip->del;   das=ip->dvi;
    _lgetdvi=(double)dvt; 
    if (!ifarg(1)) return _lgetdvi; 
    iarg=1;
    if (hoc_is_double_arg(iarg)) {
      av1=2;
      flag=*getarg(iarg++);
      getactive=(int)flag;
      flag-=(double)getactive; 
      if (flag!=0) flag=floor(flag*10+hoc_epsilon); 
    } else av1=1; 
    
    voi=vector_arg(iarg++); 
    if (flag==2) { x1=vector_newsize((IvocVect*)voi,CTYPi); for (i=0;i<CTYPi;i++) x1[i]=0;
    } else x1=vector_newsize((IvocVect*)voi,dvt);
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x2=vector_newsize((IvocVect*)voi,dvt);  a2=1; }
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x3=vector_newsize((IvocVect*)voi,dvt); a3=1;}
    if (ifarg(iarg)) { 
      voi=vector_arg(iarg++); x4=vector_newsize((IvocVect*)voi,dvt); a4=1;
      voi=vector_arg(iarg++); x5=vector_newsize((IvocVect*)voi,dvt);
    }
    idty=(double)(FOFFSET+ip->id)+1e-2*(double)ip->type+1e-3*(double)ip->inhib+1e-4;
    prty=ip->type; sy=ip->inhib?GA:AM;
    for (i=0,j=0;i<dvt;i++) {
      qp = id0ptr(das[i]->_prop); 
      if (getactive && (qp->dead || ip->sprob[i]==0)) continue;
      if (flag==1) { x1[j]=(double)qp->type; 
      } else if (flag==2) { x1[qp->type]++; 
      } else if (flag==3) { x1[j]=(double)qp->col; 
      } else x1[j]=(double)qp->id;
      if (a2) x2[j]=dbs[i];
      if (a3) x3[j]=(double)ip->sprob[i];
      if (a4) {
        poty = qp->type;
        if (seadsetting==2) { 
          for(ii=0;ii<2;ii++) y[ii]=WMAT(prty,poty,sy+ii)*WD0(prty,poty,sy+ii);
        } else {
          if (seadsetting==1) { 
            sead=(unsigned int)(FOFFSET+ip->id)*qp->id*seedstep; 
          } else { 
            hsh[0]=(double)(FOFFSET+ip->id); hsh[1]=(double)(qp->id); hsh[2]=seedstep;
            sead=hashseed2(3, hsh);
          }
          mcell_ran4(&sead, y, 2, 1.);
          for(ii=0;ii<2;ii++) {
            y[ii]=2*WVAR*(y[ii]+0.5/WVAR-0.5)*WMAT(prty,poty,sy+ii)*WD0(prty,poty,sy+ii); }
        }
        x4[j]=y[0]; x5[j]=y[1];
      }
      j++;
    }
    if (flag!=2 && j!=dvt) for (i=av1;i<iarg;i++) vector_resize(vector_arg(i),j);
    _lgetdvi=(double)j; 
  }
  ENDVERBATIM
}



FUNCTION getconv () {
VERBATIM 
{
  int iarg,i,j,k,dvt,sz,prfl,getactive; double *x,flag;
  IvocVect* voi; Point_process **das; id0 *pp;
  ip=IDP; pg=ip->pg; 
  sz=ip->dvt; 
  getactive=0;
  if (ifarg(iarg=1) && hoc_is_double_arg(iarg)) {
    flag=*getarg(iarg++);
    getactive=(int)flag;
    flag-=(double)getactive; 
    if (flag!=0) flag=floor(flag*10+hoc_epsilon);
  }
  if (!ifarg(iarg)) prfl=0; else { prfl=1;
    voi=vector_arg(iarg); 
    if (flag==2.) { x=vector_newsize(voi,CTYPi); for (i=0;i<CTYPi;i++) x[i]=0;
    } else x=vector_newsize(voi,sz); 
  } 
  for (i=0,k=0; i<cesz; i++) {
    lop(ce,i);
    if (getactive && qp->dead) continue;
    dvt=qp->dvt; das=qp->dvi;
    for (j=0;j<dvt;j++) {
      if (getactive && qp->sprob[j]==0) continue;
      if (ip == id0ptr(das[j]->_prop)) {
        if (prfl) {
          if (flag!=2.0 && k>=sz) x=vector_newsize(voi,sz*=2);
          if (flag==1.0) { x[k]=(double)qp->type; 
          } else if (flag==2.0) { x[qp->type]++; 
          } else x[k]=(double)qp->id;
        } 
        k++;
        break;
      }
    }
  }
  if (prfl && flag!=2) vector_resize(voi,k);
  _lgetconv=(double)k;
}
ENDVERBATIM
}






FUNCTION adjlist () {
  VERBATIM
  Object* pList = *hoc_objgetarg(1);
  int iListSz=ivoc_list_count(pList),iCell,iStartID=ifarg(2)?*getarg(2):0,\
    iEndID=ifarg(3)?*getarg(3):cesz-1;
  int skipinhib = ifarg(4)?*getarg(4):0, i,j,nv,*pused=(int*)calloc(cesz,sizeof(int)),iSyns=0;
  double **vvo = (double**)malloc(sizeof(double*)*iListSz),\
    *psyns=(double*)calloc(cesz,sizeof(double));
  id0* rp;
  for(iCell=iStartID;iCell<=iEndID;iCell++){
    if(verbose && iCell%1000==0) printf("%d ",iCell);
    lop(ce,iCell);
    if(!qp->dvt || (skipinhib && qp->inhib)){
      list_vector_resize(pList,iCell,0);
      continue;
    }
    iSyns=0;
    for(j=0;j<qp->dvt;j++){      
      rp = id0ptr(qp->dvi[j]->_prop); 
      if(skipinhib && rp->inhib) continue; 
      if(!rp->dead && qp->sprob[j]>0. && !pused[rp->id]){      
        pused[rp->id]=1;
        psyns[iSyns++]=rp->id;
      }
    }
    list_vector_resize(pList, iCell, iSyns);
    list_vector_px(pList, iCell, &vvo[iCell]);
    memcpy(vvo[iCell],psyns,sizeof(double)*iSyns);
    for(j=0;j<iSyns;j++)pused[(int)psyns[j]]=0;
  }
  free(vvo);  free(pused);  free(psyns);
  if (verbose) printf("\n");
  return 1.0;
  ENDVERBATIM
}

FUNCTION rddvi () {
  VERBATIM
  Point_process *pnnt;
  FILE* fp;
  int i, iCell;
  unsigned int iOutID;
  Object* lb;
  fp=hoc_obj_file_arg(1);
  printf("reading: ");
  for(iCell=0;iCell<cesz;iCell++){
    if(iCell%1000==0)printf("%d ",iCell);
    lop(ce,iCell);
    fread(&qp->id,sizeof(unsigned int),1,fp); 
    fread(&qp->type,sizeof(unsigned char),1,fp); 
    fread(&qp->col,sizeof(unsigned int),1,fp); 
    fread(&qp->dead,sizeof(unsigned char),1,fp); 
    fread(&qp->dvt,sizeof(unsigned int),1,fp); 
    
    if(qp->del){ free(qp->del); free(qp->dvi); free(qp->sprob);
      qp->dvt=0; qp->dvi=(Point_process**)0x0; qp->del=(double*)0x0; qp->sprob=(unsigned char *)0x0; }
    
    if(!qp->dvt) continue;
    qp->dvi = (Point_process**)malloc(sizeof(Point_process*)*qp->dvt);  
    for(i=0;i<qp->dvt;i++){
      fread(&iOutID,sizeof(unsigned int),1,fp); 
      if (!(lb=ivoc_list_item(ce,iOutID))) {
        printf("INTF:callback %d exceeds %d for list ce\n",iOutID,cesz); hxe(); }
      qp->dvi[i]=(Point_process *)lb->u.this_pointer;
    }
    qp->del = (double*)malloc(sizeof(double)*qp->dvt);
    fread(qp->del,sizeof(double),qp->dvt,fp); 
    qp->sprob = (unsigned char*)malloc(sizeof(unsigned char)*qp->dvt);
    fread(qp->sprob,sizeof(unsigned char),qp->dvt,fp); 
  }
  printf("\n");
  return 1.0;
  ENDVERBATIM
}

FUNCTION svdvi () {
  VERBATIM
  Point_process *pnnt;
  FILE* fp;
  int i , iCell;
  fp=hoc_obj_file_arg(1);
  printf("writing: ");
  for(iCell=0;iCell<cesz;iCell++){
    if(iCell%1000==0)printf("%d ",iCell);
    lop(ce,iCell);
    fwrite(&qp->id,sizeof(unsigned int),1,fp); 
    fwrite(&qp->type,sizeof(unsigned char),1,fp); 
    fwrite(&qp->col,sizeof(unsigned int),1,fp); 
    fwrite(&qp->dead,sizeof(unsigned char),1,fp); 
    fwrite(&qp->dvt,sizeof(unsigned int),1,fp); 
    if(!qp->dvt)continue; 
    for(i=0;i<qp->dvt;i++){
      pnnt=qp->dvi[i];
      fwrite(&(id0ptr(pnnt->_prop)->id), sizeof(unsigned int), 1, fp); 
    }
    fwrite(qp->del,sizeof(double),qp->dvt,fp); 
    fwrite(qp->sprob,sizeof(unsigned char),qp->dvt,fp); 
  }
  printf("\n"); 
  return 1.0;
  ENDVERBATIM
}









FUNCTION setdvir () {
  VERBATIM
  ListVec* pListWires,*pListDels;
  int i,dn,flag,dvt,idvfl,iCell,iStartID,iEndID,nidv,end; 
  double *y, *d, *idvec; unsigned char pdead;
  pListWires = AllocListVec(*hoc_objgetarg(1));
  idvfl=flag=0; iStartID=0; iEndID=cesz-1;
  if(!pListWires){printf("setalldvi ERRA: problem initializing wires list arg!\n"); hxe();}
  pListDels = AllocListVec(*hoc_objgetarg(2));
  if(!pListDels){ printf("setalldvi ERRA: problem initializing delays list arg!\n");
    FreeListVec(&pListWires); hxe(); }
  if (ifarg(3) && !ifarg(4)) { 
    flag=(int)*getarg(3); 
  } else if (hoc_is_double_arg(3)) {
    iStartID=(int)*getarg(3);
    iEndID = (int)*getarg(4);
    if(ifarg(5)) flag=(int)*getarg(5);
  } else {
    nidv=vector_arg_px(3, &idvec);
    idvfl=1;
    if (ifarg(4)) flag=(int)*getarg(4);
  }
  end=idvfl?nidv:(iEndID-iStartID+1);
  for (i=0;i<end;i++) {
    if(i%1000==0) printf("%d",i/1000);
    iCell=idvfl?idvec[i]:(iStartID+i);
    lop(ce,iCell);
    if (qp->dead) continue;
    y=pListWires->pv[i]; dvt=pListWires->plen[i];
    if(!dvt) continue; 
    d=pListDels->pv[i];  dn=pListDels->plen[i];
    if (dn!=dvt) {printf("setdvir() ERR vec sizes for wire,delay list entries not equal %d: %d %d\n",i,dvt,dn); hxe();}
    setdvi2(y,d,dvt,flag);
  }
  FreeListVec(&pListWires);
  FreeListVec(&pListDels);
  return 1.0;
  ENDVERBATIM
}

PROCEDURE clrdvi () {
  VERBATIM
  int i;
  for (i=0;i<cesz;i++) { 
    lop(ce,i);
    if (qp->dvt!=0x0) {
      free(qp->dvi); free(qp->del); free(qp->sprob);
      qp->dvt=0; qp->dvi=(Point_process**)0x0; qp->del=(double*)0x0; qp->sprob=(unsigned char *)0x0;
    }
  }
  ENDVERBATIM
}


FUNCTION setdviv () {
  VERBATIM
  int i,j,k,nprv,dvt; double *prv,*pov,*dlv,x;
  nprv=vector_arg_px(1, &prv);
  i=vector_arg_px(2, &pov);
  j=vector_arg_px(3, &dlv);
  if (nprv!=i || i!=j) {printf("intf:setdviv ERRA: %d %d %d\n",nprv,i,j); hxe();}
  
  if (scrsz<cesz) scrset(cesz); 
  for (i=0;i<cesz;i++) scr[i]=0;
  for (i=0,j=-1;i<nprv;i++) {
    if (j>(int)prv[i]){printf("intf:setdviv ERRA vecs should be sorted by prid vec\n");hxe();}
    j=(int)prv[i];
    scr[j]++;
  }
  for (i=0,x=-1,k=0;i<nprv;i+=dvt) { if(i%1000==0) printf(".");
    if (prv[i]!=x) lop(ce,(unsigned int)(x=prv[i]));
    if (qp->dead) continue;
    dvt=scr[(int)x]; 
    setdvi2(pov+k,dlv+k,dvt,1);
    k+=dvt;
  }
  return (double)k;
  ENDVERBATIM
}


VERBATIM
static void finishdvi2 (struct ID0* p) {
  Point_process **da,**das;
  double *db,*dbs;
  int i, dvt;
  db=p->del;
  da=p->dvi; 
  dvt=p->dvt;
  dbs=(double*)malloc(dvt*sizeof(double)); 
  das=(Point_process**)malloc(dvt*sizeof(Point_process*)); 
  gsort2(db,da,dvt,dbs,das);
  p->del=dbs; p->dvi=das; 
  free(db); free(da);
  p->sprob=(unsigned char*)realloc((void*)p->sprob,(size_t)dvt*sizeof(char));
  for (i=0;i<dvt;i++) p->sprob[i]=1; 
}
ENDVERBATIM


PROCEDURE finishdvir () {
  VERBATIM
  int iCell;
  for(iCell=0;iCell<cesz;iCell++){
    lop(ce,iCell);
    finishdvi2(qp);
  }
  ENDVERBATIM
}


PROCEDURE finishdvi () {
VERBATIM
  finishdvi2(IDP);
ENDVERBATIM
}


PROCEDURE setdvi () {
VERBATIM {
  int i,dvt,flag; double *d, *y;
  if (! ifarg(1)) {printf("setdvi(v1,v2[,flag]): v1:cell#s; v2:delays\n"); return 0; }
  ip=IDP; pg=ip->pg; 
  if (ip->dead) return 0;
  dvt=vector_arg_px(1, &y);
  i=vector_arg_px(2, &d);
  if (ifarg(3)) flag=(int)*getarg(3); else flag=0;
  if (i!=dvt || i==0) {printf("setdvi() ERR vec sizes: %d %d\n",dvt,i); hxe();}
  setdvi2(y,d,dvt,flag);
  }
ENDVERBATIM
}

VERBATIM


static void setdvi2 (double *y,double *d,int dvt,int flag) {
  int i,j,ddvi; double *db, *dbs; unsigned char pdead; unsigned int b,e;
  Object *lb; Point_process *pnnt, **da, **das;
  ddvi=(int)DEAD_DIV;
  ip=IDP; 
  if (flag==0) { b=0; e=dvt; 
    if (ip->dvi) { 
      free(ip->dvi); free(ip->del); free(ip->sprob); 
      ip->dvt=0; ip->dvi=(Point_process**)0x0; ip->del=(double*)0x0; ip->sprob=(unsigned char *)0x0;
    } 
  } else { 
    if (ip->dvt==0) {ip->dvi=(Point_process**)0x0; ip->del=(double*)0x0; ip->sprob=(unsigned char *)0x0;}
    b=ip->dvt; 
    e=ip->dvt+dvt; 
  }
  da=(Point_process **)realloc((void*)ip->dvi,(size_t)(e*sizeof(Point_process *)));
  db=(double*)realloc((void*)ip->del,(size_t)(e*sizeof(double)));
  for (i=b,j=0;j<dvt;j++) { 
    
    if (!(lb=ivoc_list_item(ce,(unsigned int)y[j]))) {
      printf("INTF:callback %g exceeds %d for list ce\n",y[j],cesz); hxe(); }
      pnnt=(Point_process *)lb->u.this_pointer;
      if (ddvi==1 || !(pdead = id0ptr(pnnt->_prop)->dead)) {
        da[i]=pnnt; db[i]=d[j]; i++;
      }
  }
  if ((dvt=i)<e) { 
    da=(Point_process **)realloc((void*)da,(size_t)(e*sizeof(Point_process *)));
    db=(double*)realloc((void*)db,(size_t)(e*sizeof(double)));
  }
  ip->dvt=dvt; ip->del=db; ip->dvi=da;
  if (flag!=1) finishdvi2(ip); 
}
ENDVERBATIM



PROCEDURE prune () {
  VERBATIM 
  {
  id0* ppost; double *x, p; int nx,j,potype;
  ip=IDP; pg=ip->pg;
  if (hoc_is_double_arg(1)) { 
    p=*getarg(1);
    if (p<0 || p>1) {printf("INTF:pruneERR0:need # [0,1] to prune [ALL,NONE]: %g\n",p); hxe();}
    if (p==1.) printf("INTFpruneWARNING: pruning 100% of cell %d\n",ip->id);
    if (verbose && ip->dvt>dscrsz) {
      printf("INTFpruneB:Div exceeds dscrsz: %d>%d\n",ip->dvt,dscrsz); hxe(); }
    if (p==0.) {
      for (j=0;j<ip->dvt;j++) ip->sprob[j]=1; 
      return 0; 
    }
    potype=ifarg(2)?(int)*getarg(2):-1;
    sead=(ifarg(3))?(unsigned int)*getarg(3):GetDVIDSeedVal(ip->id);
    mcell_ran4(&sead, dscr , ip->dvt, 1.0); 
    if(potype==-1){ 
      for (j=0;j<ip->dvt;j++) if (dscr[j]<p) ip->sprob[j]=0; 
    } else { 
      for (j=0;j<ip->dvt;j++){
        ppost = id0ptr(ip->dvi[j]->_prop); 
        if (ppost->type==potype && dscr[j]<p) ip->sprob[j]=0; 
      }
    }
  } else { 
    if (verbose) printf("INTF WARNING prune(vec) deprecated: use intf.sprob(vec) instead\n");
    nx=vector_arg_px(1,&x);
    if (nx!=ip->dvt) {printf("INTF:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=(unsigned char)x[j];
  }
  }
ENDVERBATIM
}

PROCEDURE sprob () {
  VERBATIM 
  {
  double *x; int nx,j;
  ip=IDP; pg=ip->pg;
  nx=vector_arg_px(1,&x);
  if (nx!=ip->dvt) {printf("INTF:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
  if (ifarg(2)) { 
    if (!hoc_is_str_arg(2)) { printf("INTF sprob()ERRA: only legit 2nd arg is 'GET'\n"); hxe();
    } else for (j=0;j<ip->dvt;j++) x[j]=(double)ip->sprob[j];
  } else {
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=(unsigned char)x[j];
  }
  }
ENDVERBATIM
}



PROCEDURE turnoff () {
  VERBATIM {
  int nx,ny,i,j,k,dvt; double poid,*x,*y; Point_process **das; unsigned char off;
  nx=vector_arg_px(1,&x);
  ny=vector_arg_px(2,&y);
  if (ifarg(3)) off=(unsigned char)*getarg(3); else off=0;
  for (i=0;i<nx;i++) { 
    lop(ce,(unsigned int)x[i]); 
    dvt=qp->dvt; das=qp->dvi;
    for (j=0;j<dvt;j++) {
      ip = id0ptr(das[j]->_prop); 
      poid=(double)ip->id; 
      for (k=0;k<ny;k++) {
        if (poid==y[k]) {
          qp->sprob[j]=off; break;
        }
      }
    }
  }
  }
  ENDVERBATIM
}

VERBATIM 

void gsort2 (double *db, Point_process **da,int dvt,double *dbs, Point_process **das) {
  int i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, (int*)scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
  }
}
ENDVERBATIM

PROCEDURE freedvi () {
  VERBATIM
  { 
    int i, poty; id0 *jp;
    jp=IDP;
    if (jp->dvi) {
      free(jp->dvi); free(jp->del); free(jp->sprob);
      jp->dvt=0; jp->dvi=(Point_process**)0x0; jp->del=(double*)0x0; jp->sprob=(unsigned char *)0x0;
    }
  }
  ENDVERBATIM
}


PROCEDURE pgset () {
  VERBATIM
  ip->pg=pg; 
  ENDVERBATIM
}

FUNCTION qstats () {
  VERBATIM {
    double stt[3]; int lct,flag; FILE* tfo;
    if (ifarg(1)) {tfo=hoc_obj_file_arg(1); flag=1;} else flag=0;
    lct=cty[IDP->type];
    _lqstats = nrn_event_queue_stats(stt);
    printf("SPIKES: %d (%ld:%ld)\n",IDP->spkcnt,spikes[lct],blockcnt[lct]);
    printf("QUEUE: Inserted %g; removed %g\n",stt[0],stt[2]);
    if (flag) {
      fprintf(tfo,"SPIKES: %d (%ld:%ld);",IDP->spkcnt,spikes[lct],blockcnt[lct]);
      fprintf(tfo,"QUEUE: Inserted %g; removed %g remaining: %g\n",stt[0],stt[2],_lqstats);
    }
  }
  ENDVERBATIM
}

FUNCTION qsz () {
  VERBATIM {
    double stt[3];
    _lqsz = nrn_event_queue_stats(stt);
  }
  ENDVERBATIM
}

PROCEDURE qclr () {
  VERBATIM {
    clear_event_queue();
  }
  ENDVERBATIM
}


PROCEDURE jitcondiv () {
  VERBATIM {
  Symbol *sym; int i,j; unsigned int sz; char *name="ce";
  pg=(postgrp *)malloc(sizeof(postgrp));
  if (ifarg(1)) name = gargstr(1);
  sym = hoc_lookup(name); ce = (*(hoc_objectdata[sym->u.oboff].pobj));
  sym = hoc_lookup("CTYP"); CTYP = (*(hoc_objectdata[sym->u.oboff].pobj));
  if (installed==2.0) { 
    sz=ivoc_list_count(ce);
    if (sz==cesz) printf("\t**** INTF WARNING cesz unchanged: INTF(s) created off-list ****\n");
  } else installed=2.0;
  cesz = ivoc_list_count(ce);
  cty[0]=DP; cty[1]=SU; cty[2]=IN; 
  CTYPi=HVAL("CTYPi"); STYPi=HVAL("STYPi"); dscrsz=HVAL("scrsz");
  pg->ix =HPTR("ix"); pg->ixe=HPTR("ixe"); 
  pg->dvg=HPTR("div"); pg->numc=HPTR("numc"); 
  pg->wmat=HPTR("wmat"); pg->wd0=HPTR("wd0");
  pg->delm=HPTR("delm"); pg->deld=HPTR("deld");
  dscr=HPTR("scr");
  if (!ce) {printf("INTF jitcondiv ERRA: ce not found\n"); hxe();}
  if (ivoc_list_count(CTYP)!=CTYPi){
    printf("INTF jitcondiv ERRB: %d %d\n",ivoc_list_count(CTYP),CTYPi); hxe(); }
  for (i=0;i<cesz;i++) { lop(ce,i); qp->pg=pg; } 
  
  printf("Checking for possible seg error in double arrays: CTYPi==%d: ",CTYPi);
  
  printf("%d %d %d ",DVG(CTYPi-1,CTYPi-1),(int)pg->ix[CTYPi-1],(int)pg->ixe[CTYPi-1]);
  printf("%g %g ",WMAT(CTYPi-1,CTYPi-1,STYPi-1),WD0(CTYPi-1,CTYPi-1,STYPi-1));
  printf("%g %g ",DELM(CTYPi-1,CTYPi-1),DELD(CTYPi-1,CTYPi-1));
  printf("%d %g\n",dscrsz,dscr[dscrsz-1]); 
  for (i=0,j=0;i<CTYPi;i++) if (ctt(i,&name)!=0) {
    cty[j]=i; CNAME[j]=name;
    j++;
    if (j>=CTYPp) {printf("jitcondiv() INTERRA\n"); hxe();}
  }
  CTYN=j; 
  for (i=0;i<CTYN;i++) printf("%s(%d)=%g ",CNAME[i],cty[i],NUMC(cty[i]));
  printf("\n%d cell types being used\n",CTYN);
  }
  ENDVERBATIM  
}


PROCEDURE jitrec () {
  VERBATIM {
  int i;
  if (! ifarg(2)) { 
    jrmax=0; jridv=0x0; jrtvv=0x0;
    return 0;
  }
  i =   vector_arg_px(1, &jrid); 
  jrmax=vector_arg_px(2, &jrtv);
  jridv=vector_arg(1); jrtvv=vector_arg(2);
  jrmax=vector_buffer_size(jridv);
  if (jrmax!=vector_buffer_size(jrtvv)) {
    printf("jitrec() ERRA: not same size: %d %d\n",i,jrmax); jrmax=0; hxe(); }
  jri=jrj=0; 
  }
  ENDVERBATIM
}


FUNCTION scsv () {
  VERBATIM {
  int ty=4; int i,j; unsigned int cnt=0;
  name = gargstr(1);
  if ( !(wf1 = fopen(name,"w"))) { printf("Can't open %s\n",name); hxe(); }
  fwrite(&cesz,sizeof(int),1,wf1);
  fwrite(&ty,sizeof(int),1,wf1);
  for (i=0,j=0;i<cesz;i++,j++) { 
    lop(ce,i); 
    if (qp->spkcnt) {
      dscr[j]=(double)(qp->spkcnt); 
      cnt++;
    } else dscr[j]=0.0;
    if (j>=dscrsz) {
      fwrite(dscr,(size_t)sizeof(double),(size_t)dscrsz,wf1);
      fflush(wf1);
      j=0;
    }
  }
  if (j>0) fwrite(dscr,(size_t)sizeof(double),(size_t)j,wf1);
  fclose(wf1);
  _lscsv=(double)cnt;
  }
  ENDVERBATIM
}



FUNCTION spkcnt () {
  VERBATIM {
  int nx, ny, i,j, ix, c, min, max, flag; unsigned int sum; double *y,*x;
  nx=ny=min=max=flag=0; i=1;
  if (ifarg(i)) {
    if (hoc_is_object_arg(i)) { 
      ny = vector_arg_px(i, &y); i++;
    } else if (ifarg(i+1)) {
      min=(int)*getarg(i); max=(int)*getarg(i+1); i+=2;
    }
  }
  while (ifarg(i)) { 
    if (hoc_is_object_arg(i)) { 
      nx = vector_arg_px(i, &x);
    } else flag=(int)*getarg(i);
    i++;
  }
  if (ny) max=ny; else if (max==0) max=cesz; else max+=1; 
  if (nx && nx!=max-min) {
    printf("INTF spkcnt() ERR: Vectors not same size %d %d\n",nx,max-min);hxe();}
  for  (i=min, sum=0;i<max;i++) { 
    if (ny) lop(ce,(int)y[i]); else lop(ce,i);
    if (flag==2) sum+=(c=qp->blkcnt); else sum+=(c=qp->spkcnt);
    if (nx) x[i]=(double)c;
    if (flag==1) qp->spkcnt=qp->blkcnt=0;
  }
  _lspkcnt=(double)sum;
  }
  ENDVERBATIM
}


PROCEDURE probejcd () {
  VERBATIM {  int i,a[4];
    for (i=1;i<=3;i++) a[i]=(int)*getarg(i);
    printf("CTYPi: %d, STYPi: %d, ",CTYPi,STYPi);
    
    printf("wmat: %g, wd0: %g\n",WMAT(a[1],a[2],a[3]),WD0(a[1],a[2],a[3]));
  }
  ENDVERBATIM  
}


PROCEDURE randspk () {
  VERBATIM 
  ip=IDP;  
  if (ip->rvi > ip->rve) { 
    ip->input=0;           
    nxt=-1.;
  } else if (t==0) {     
    nxt=vsp[ip->rvi];
    WEX=wsp[ip->rvi++];
  } else {     
    while ((nxt=vsp[ip->rvi++]-t)<=1e-6) { 
      if (ip->rvi-1 > ip->rve) { printf("randspk() ERRA: "); chk(2.); hxe(); }
    }
    WEX=wsp[ip->rvi-1]; 
  }
  ENDVERBATIM
  
}


PROCEDURE vers () {
  printf("$Id
}


VERBATIM
void val (double xx, double ta) { 
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
  vii[7]=VTH + (VTHR-VTH)*EXP(-(xx-trrs)/tauRR);
}
ENDVERBATIM


VERBATIM
void valps (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  vii[4]=esinr(xx-tGB);
  
  vii[6]=vii[1]+vii[2]-vii[3];
}
ENDVERBATIM


PROCEDURE record () {
  VERBATIM {
  int i,j,k,nz; double ti;
  vp = SOP;
  if (tg>=t) return 0;
  if (ip->record==1) {
    while ((int)vp->p >= (int)vp->size-(int)((t-tg)/vdt)-10) { 
      vp->size*=2;
      for (k=0;k<NSV;k++) if (vp->vv[k]!=0x0) vp->vvo[k]=vector_newsize(vp->vv[k], vp->size);
      
    }
  } else if ((int)vp->p > (int)vp->size-(int)((t-tg)/vdt)) { 
    nz=(int)((t-tg)/vdt);
    for (k=0;k<NSV;k++) if (vp->vv[k]!=0x0) {
      if (nz>vp->size) {pid(); printf("Record WARNING: vec too short: %d %d\n",nz,vp->size);
        vp->p=0;
      } else {
        for (i=nz,j=0; i<vp->size; i++,j++) vp->vvo[k][j]=vp->vvo[k][i];
        vp->p=vp->size-nz;
      }
    }
  }
  for (ti=tg;ti<=t && vp->p < vp->size;ti+=vdt,vp->p++) { 
    val(ti,tg);  
    if (vp->vvo[0]!=0x0) vp->vvo[0][vp->p]=ti;
    for (k=1;k<NSV-1;k++) if (vp->vvo[k]!=0x0) { 
      vp->vvo[k][vp->p]=vii[k]+RMP;
    }
    for (;k<NSV;k++) if (vp->vvo[k]!=0x0) { 
      vp->vvo[k][vp->p]=vii[k]; 
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
  if (vp->p > vp->size || vp->vvo[6]==0) return 0;
  if (vp->vvo[0]!=0x0) vp->vvo[0][vp->p-1]=_lx;
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
    double last,lstt;
    if (! ifarg(1)) {printf("Return initvspks(ivspks,vspks,wvspks)\n"); return 0.;}
    if (isp!=NULL) clrvspks();
    ip=IDP;  err=0;
    i = vector_arg_px(1, &isp); 
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("initvspks ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("initvspks ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("initvspks ERR: 3rd vec is of different size\n");}
    vspn=max;
    if (!ce) {printf("Need global ce for initvspks() since intf.mod501\n"); hxe();}
    for (i=0,last=-1; i<max; ) { 
      if (isp[i]!=last) { 
        lop(ce,(unsigned int)isp[i]);
        qp->rvb=qp->rvi=i;
        qp->vinflg=1;
        last=isp[i];
        lstt=vsp[i];
        i++;
      }
      for (; i<max && isp[i] == last; i++) { 
        if (vsp[i]<=lstt) { err=1; 
          printf("initvspks ERR: nonmonotonic for cell#%d: %g %g\n",qp->id,lstt,vsp[i]); }
          lstt=vsp[i];
      }
      qp->rve=i-1;
      if (subsvint>0) { 
        vsp[qp->rve] = vsp[qp->rvb]+subsvint;
        wsp[qp->rve] = wsp[qp->rvb];
      }
      if (err) { qp->rve=0; hxe(); }
    }
  }
  ENDVERBATIM
}




PROCEDURE shock () {
  VERBATIM 
  {int max, i,err;
    double last, lstt, *isp, *vsp, *wsp;
    if (! ifarg(1)) {printf("Return shock(ivspks,vspks,wvspks)\n"); return 0.;}
    ip=IDP;  err=0;
    i = vector_arg_px(1, &isp); 
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("shock ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("shock ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("shock ERR: 3rd vec is of different size\n");}
    vspn=max;
    if (!ce) {printf("Need global ce for shock()\n"); hxe();}
    for (i=0,last=-1; i<max; ) { 
      if (isp[i]!=last) { 
        lop(ce,(unsigned int)isp[i]);
        WEX=-1e9; 
  #if defined(t)
        net_send((void**)0x0, wts,pmt,t+vsp[i],2.0); 
  #else
        net_send((void**)0x0, wts,pmt,vsp[i],2.0); 
  #endif
        i++;
      }
    }
  }
  ENDVERBATIM
}

PROCEDURE clrvspks () {
 VERBATIM {
 unsigned int i;
 for (i=0; i<cesz; i++) {
   lop(ce,i);
   qp->vinflg=0;
 }   
 }
 ENDVERBATIM
}




PROCEDURE trvsp ()
{
  VERBATIM 
  int i, flag; 
  double ind, local_t0;
  ip=IDP;
  flag=(int) *getarg(1);
  if (subsvint==0.) {printf("trvsp"); return(0.);}
  ind=isp[0];
  local_t0=vsp[0];
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
        vsp[i-1] = local_t0 + subsvint;
        ind = isp[i];
        local_t0 = vsp[i];
      }
    }
    vsp[vspn-1] = local_t0 + subsvint;
  } else {printf("trvsp flag %d not recognized\n",flag); hxe();}
  ENDVERBATIM
}




PROCEDURE initjttr () {
  VERBATIM 
  {int max, i, err=0;
    jtpt=0;
    if (! ifarg(1)) {printf("Return initjttr(vec)\n"); return(0.);}
    max=vector_arg_px(1, &jsp);
    if (max==0) {err=1; printf("initjttr ERR: vec not initialized\n");}
    for (i=0; i<max; i++) if (jsp[i]<=0) {err=1;
      printf("initjttr ERR: vec should be >0: %g\n",jsp[i]);}
    if (err) { jsp=nil; jtmax=0.; return(0.); }
    if (max != jtmax) {
      printf("WARNING: resetting jtmax_INTF to %d\n",max); jtmax=max; }
  }
  ENDVERBATIM
}


VERBATIM


static void lop (Object *ob, unsigned int i) {
  Object *lb;
  lb = ivoc_list_item(ob, i);
  if (! lb) { printf("INTF:lop %d exceeds %d for list ce\n",i,cesz); hxe();}
  pmt=ob2pntproc(lb);
  qp = id0ptr(pmt->_prop); 
  
#ifdef NRN_MECHANISM_DATA_IS_SOA
  neuron::legacy::set_globals_from_prop(pmt->_prop, _ml_real, _ml, _iml);
#else
  _p = pmt->_prop->param;
#endif
  _ppvar = get_dparam(pmt->_prop);
}


void stoppo () {
}


static int ctt (unsigned int i, char** name) {
  Object *lb;
  if (NUMC(i)==0) return 0; 
  lb = ivoc_list_item(CTYP, i);
  if (! lb) { printf("INTF:ctt %d exceeds %d for list CTYP\n",i,CTYPi); hxe();}
  {*name=*(lb->u.dataspace->ppstr);}
  return (int)NUMC(i);
}
ENDVERBATIM


PROCEDURE test () {
  VERBATIM
  char *str; int x;
  x=ctt(7,&str); 
  printf("%s (%d)\n",str,x);
  ENDVERBATIM
}


PROCEDURE lof () {
VERBATIM {
  Object *ob; int num,i,ii,j,k,si,nx;  double *vvo[7], *par; IvocVect *vv[7];
  ob = *(hoc_objgetarg(1));
  si=(int)*getarg(2);
  num = ivoc_list_count(ob);
  if (num!=7) { printf("INTF lof ERR %d>7\n",num); hxe(); }
  for (i=0;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);
    if (i==0) nx=j;
    if (j!=nx) { printf("INTF lof ERR %d %d\n",j,nx); hxe(); }
  }
  
  
  
  
  
  
  
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
    printf("INTF record ERR: pointer not initialized\n"); hoc_execerror("",0); 
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
  printf("ID:%d; typ: %d; rec:%d wrec:%d inp:%d jtt:%d invl:%d\n",ip->id,ip->type,ip->record,ip->wrec,ip->input,ip->jttr,ip->invl0);
  if (lfg==1) {
    if (SOP!=nil) {
      vp = SOP;
      printf("p %d size %d tg %g\n",vp->p,vp->size,tg);
      for (i=0;i<NSV;i++) if (vp->vv[i]) printf("%d %x %x;",i,vp->vv[i],vp->vvo[i]);
    } else printf("Recording pointers not initialized");
  }
  if (lfg==2) { 
    printf("Global vectors for input and jitter (jttr): \n");
    if (vsp!=nil) printf("VSP: %x (%d/%d-%d)\n",vsp,ip->rvi,ip->rvb,ip->rve); else printf("no VSP\n");
    if (jsp!=nil) printf("JSP: %x (%d/%d)\n",jsp,jtpt,jtmax); else printf("no JSP\n");
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
    for (i=0;i<NSW;i++) printf("%d %p %p;",i,ww[i],wwo[i]);
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
  if (ifarg(1)) IDP->id = (unsigned int) *getarg(1);
  _lid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION type () {
  VERBATIM
  if (ifarg(1)) IDP->type = (unsigned char) *getarg(1);
  _ltype = (double)IDP->type;
  ENDVERBATIM
}

FUNCTION col () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->col = (unsigned int) *getarg(1);
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
  {int i;
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
  } else if (strcmp(name,"VTHC")==0) { _li=7.;
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
  else if (i==7) printf("VTHC\n");
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
  {int i,j,k,max,wrp; double ti,scale;
  for (i=0;i<WRNUM && (wrp=(int)IDP->wreci[i])>-1;i++) {
    
    scale=(double)IDP->wscale[i];
    if (_lte<1.e9) { 
      if (scale>0) {
        max=-(int)_tmin_popspk; 
        k=(int)floor((_lte-rebeg)/vdt+0.5);
        for (j= -max;j<=max && k+j>0 && k+j<wwsz;j++) {
          wwo[wrp][k+j] += scale*_t_Psk[j+max]; 
        }
      }
    } else if (twg>=t) { return 0;
    } else {
      for (ti=twg,k=(int)floor((twg-rebeg)/vdt+0.5);ti<=t && k<wwsz;ti+=vdt,k++) { 
        valps(ti,twg);  
        wwo[wrp][k]+=vii[6];
        if (IDP->dbx==-1) printf("%g:%g ",vii[6],wwo[wrp][k]);
      }
    }
  }
  if (_lte==1.e9) twg=ti;
  }
  ENDVERBATIM
}





FUNCTION wrec () {
  VERBATIM
  { int k,ix;
  ip=IDP; 
  if (ifarg(1)) {
    ix=(int)*getarg(1);
    if (ix>=1) {
      if (ix-1>=nsw) {
        printf("Attempt to save into ww[%d] but only have %d\n",ix-1,(int)nsw); hxe();}
      ip->wrec=1;
      ip->wreci[0]=(char)ix-1;
      ip->wscale[0]=1.; 
      if (ifarg(2)) ip->wscale[0]= (float)*getarg(2); 
    } else if (ix<=0) {
      ip->wrec=0;
      for (k=0;k<WRNUM;k++) { ip->wreci[k]=-1; ip->wscale[k]=-1.0; }
    } else {printf("INTF wrec ERR flag(0/1) %d\n",ip->wrec); hxe();
    }
  }
  _lwrec=(double)ip->wrec;
  }
  ENDVERBATIM
}





FUNCTION wrc () {
  VERBATIM
  { int i,ix;
  ip=IDP; 
  if (ifarg(1)) {  
    ix=(int)*getarg(1);
    if (ix<0) {
      ip->wrec=0;
      for (i=0;i<WRNUM;i++) { ip->wreci[i]=-1; ip->wscale[i]=-1.0; }
    } else {
      for (i=0;i<WRNUM && ip->wreci[i]!=-1 && ip->wreci[i]!=ix;i++) {};
      if (i==WRNUM) {
        pid(); printf("INFT wrc() ERR: out of wreci pointers (max %d)\n",WRNUM); hxe();}
      if (ix>=nsw) {printf("Attempt to save into ww[%d] but only have %d\n",ix,(int)nsw); hxe();}
      ip->wrec=1; 
      ip->wreci[i]=ix;
      if (ifarg(2)) ip->wscale[i]=(float)*getarg(2); else ip->wscale[i]=1.0;
    }
  } else {
    for (i=0;i<WRNUM;i++) printf("%d:%g ",ip->wreci[i],ip->wscale[i]);
    printf("\n");
  }
  _lwrc=(double)ip->wrec;
  }
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


FUNCTION jttr () {
  VERBATIM 
  if (jtmax>0 && jtpt>=jtmax) {  
    jtpt=0;
    printf("Warning, cycling through jttr vector at t=%g\n",t);
  }
  if (jtmax>0) _ljttr = jsp[jtpt++]; else _ljttr=0;
  ENDVERBATIM
}


PROCEDURE global_init () {
  popspk(0) 
  VERBATIM 
  { int i,j,k; double stt[3];
  if (jridv) { jri=jrj=0; vector_resize(jridv, jrmax); vector_resize(jrtvv, jrmax); }
  if (nsw>0. && wwo[0]!=0) { 
    printf("Initializing ww to record for %g (%g)\n",vdt*wwsz,vdt);
    wwpt=0;
    for (k=0;k<(int)nsw;k++) {
      vector_resize(ww[k], wwsz);
      for (j=0;j<wwsz;j++) wwo[k][j]=0.;
    }
  }
  spktot=0;
  jtpt=0;
  eventtot=0;
  errflag=0;
  for (i=0;i<CTYN;i++) blockcnt[cty[i]]=spikes[cty[i]]=0;
  }
  ENDVERBATIM
}

PROCEDURE global_fini () {
  VERBATIM
  int k;
  for (k=0;k<(int)nsw;k++) vector_resize(ww[k], (int)floor(t/vdt+0.5));
  if (jridv && jrj<jrmax) {
    vector_resize(jridv, jrj); 
    vector_resize(jrtvv, jrj);
  }
  ENDVERBATIM
}


FUNCTION fflag () { fflag=1 }
FUNCTION thrh () { thrh=VTH-RMP }

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




FUNCTION flag () {
  VERBATIM
  char *sf; static int ix,fi,setfl,nx; static unsigned char val; static double *x, delt;
  if (FLAG==OK) { 
    FLAG=0.;
    if (stoprun) {slowset=0; return 0;}
    if (IDP->dbx==-1)printf("slowset fi:%d ix:%d ss:%g delt:%g t:%g\n",fi,ix,slowset,delt,t);
    if (t>slowset || ix>=cesz) {  
      printf("Slow-setting of flag %d finished at %g: (%d,%g,%g)\n",fi,t,ix,delt,slowset); 
      slowset=0.; return 0;
    }
    if (ix<cesz) {
      lop(ce,ix);
      (&qp->type)[fi]=((fi>=iflneg)?(char)x[ix]:(unsigned char)x[ix]);
      ix++;
      #if defined(t)
      net_send((void**)0x0, wts,tpnt,t+delt,OK); 
      #else
      net_send((void**)0x0, wts,tpnt,delt,OK);
      #endif
    }
    return 0;
  }  
  if (slowset>0 && ifarg(3)) {
    printf("INTF flag() slowset ERR; attempted set during slowset: fi:%d ix:%d ss:%g delt:%g t:%g",\
           fi,ix,slowset,delt,t); 
    return 0;
  }
  ip = IDP; setfl=ifarg(3); 
  if (ifarg(4)) { slowset=*getarg(4); delt=slowset/cesz; slowset+=t; } 
  sf = gargstr(1);
  for (fi=0;fi<iflnum && strncmp(sf, &iflags[fi*4], 3)!=0;fi++) ; 
  if (fi==iflnum) {printf("INTF ERR: %s not found as a flag (%s)\n",sf,iflags); hxe();}
  if (ifarg(2)) {
    if (hoc_is_double_arg(2)) {  
      val=(unsigned char)*getarg(2);
      if (slowset) { 
        printf("NOT IMPLEMENTED\n"); 
      } else if (setfl) { 
        for (ix=0;ix<cesz;ix++) { lop(ce,ix); (&qp->type)[fi]=val; }
      } else {  
        (&ip->type)[fi]=((fi>=iflneg)?(char)val:val);
      }
    } else {
      nx=vector_arg_px(2,&x);
      if (nx!=cesz) {
        if (setfl) { printf("INTF flag ERR: vec sz mismatch: %d %d\n",nx,cesz); hxe();
        } else       x=vector_newsize(vector_arg(2),cesz);
      }
      if (setfl && slowset) { 
        ix=0;
        lop(ce,ix);
        (&qp->type)[fi]=((fi>=iflneg)?(char)x[ix]:(unsigned char)x[ix]);
        ix++;
        #if defined(t)
        net_send((void**)0x0, wts,tpnt,t+delt,OK); 
        #else
        net_send((void**)0x0, wts,tpnt,delt,OK);
        #endif
      } else for (ix=0;ix<cesz;ix++) { 
        lop(ce,ix); 
        if (setfl) { (&qp->type)[fi]=((fi>=iflneg)?(char)x[ix]:(unsigned char)x[ix]);
        } else {
          x[ix]=(double)((fi>=iflneg)?(char)(&qp->type)[fi]:(unsigned char)(&qp->type)[fi]);
        }
      }
    }
  }
  _lflag=(double)((fi>=iflneg)?(char)(&ip->type)[fi]:(unsigned char)(&ip->type)[fi]);
  ENDVERBATIM
}

FUNCTION allspck () {
  VERBATIM
  int i; double *x, sum; IvocVect* voi;
  ip = IDP;
  voi=vector_arg(1);  x=vector_newsize(voi,cesz);
  for (i=0,sum=0;i<cesz;i++) { lop(ce,i); 
    x[i]=spck;
    sum+=spck;
  }
  _lallspck=sum;
  ENDVERBATIM
}


PROCEDURE resetall () {
  VERBATIM
  int ii,i; unsigned char val;
  for (i=0;i<cesz;i++) { 
    lop(ce,i);
    Vm=RMP; VAM=0; VNM=0; VGA=0; VGB=0; VGBa=0; offsetGB=0; AHP=0; rebob=-1e9; invlt=-1; 
    t0=t; tGB=t; trrs=t; twg = t; cbur=0; spck=0; refractory=0; VTHC=VTHR=VTH; 
  }
  ENDVERBATIM
}


FUNCTION floc () {
  VERBATIM
  double x,y,z,r,min,rad, *ix; int ii,i,n,cnt; IvocVect* voi;
  cnt=0; n=1000; r=-1;
  ip = IDP;
  x = *getarg(1);
  y = *getarg(2);
  z= ifarg(3)?(*getarg(3)):1e9;
  if (ifarg(5)) {
    r= *getarg(4);
    voi=vector_arg(5);
    ix=vector_newsize(voi,n);
  } 
  for (i=0,min=1e9,ii=-1;i<cesz;i++) { 
    lop(ce,i); 
    rad=(x-xloc)*(x-xloc)+(y-yloc)*(y-yloc)+(z==1e9?0.:((z-zloc)*(z-zloc))); 
    if (r>0 && rad<r*r) {
      if (cnt>=n) ix=vector_newsize(voi,n*=2);
      ix[cnt]=(double)i;
      cnt++;
    }
    if (rad<min) { min=rad; ii=i; }
  }
  if (r>0) ix=vector_newsize(voi,cnt);
  _lfloc=(double)ii;
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

PROCEDURE rates (vv) {
  TABLE Bb DEPEND mg FROM -100 TO 50 WITH 300
  
  Bb = 1 / (1 + exp(0.062 (/mV) * -vv) * (mg / 3.57 (mM)))
}

PROCEDURE coop (x) {
  TABLE Gn DEPEND GPkd FROM 0 TO 10 WITH 100
  
  Gn = (x^4)/(x^4+GPkd) 
}