NEURON {
  ARTIFICIAL_CELL INTF6
  RANGE VAM, VNM, VGA, AHP           
  RANGE VAM2, VNM2, VGA2                  
  RANGE Vm                                
  
  RANGE tauAM, tauNM, tauGA            
  RANGE tauAM2, tauNM2, tauGA2         
  RANGE tauahp, ahpwt                  
  RANGE tauRR , RRWght                 
  RANGE RMP,VTH,Vblock,VTHC,VTHR       
  RANGE incRR 
  RANGE nbur,tbur,refrac,AHP2REF        
  RANGE invl,oinvl,WINV,invlt           
  RANGE Vbrefrac                        
  RANGE STDAM, STDNM, STDGA             
                                        
                                        
                                        
                                        
  RANGE mg0                             
  RANGE maxnmc                          
  GLOBAL EAM, ENM, EGA,mg               
  GLOBAL spkht, wwwid,wwht              
  GLOBAL stopoq                         
  
  POINTER sop                          
  RANGE  spck,xloc,yloc,zloc
  RANGE  t0,tg,twg,refractory,trrs 
  RANGE  cbur                         
  RANGE  WEX                          
  RANGE  EXSY                         
  RANGE  lfpscale                     
  GLOBAL vdt,nxt,RES,ESIN,Psk      
  GLOBAL prnum, nsw, rebeg             
  GLOBAL subsvint, jrsvn, jrsvd, jrtime, jrtm 
  GLOBAL DEAD_DIV, seedstep            
  GLOBAL seaddvioff                    
  GLOBAL WVAR,DELMIN
  GLOBAL savclock,slowset,FLAG  
  GLOBAL tmax,installed,verbose        
  GLOBAL pathbeg,pathend,PATHMEASURE,pathidtarg,pathtytarg,seadsetting,pathlen
  GLOBAL maxplastt 
  GLOBAL plaststartT 
  GLOBAL plastendT   
  GLOBAL resetplast  
  GLOBAL wsetting 
  GLOBAL ESTDP 
  GLOBAL ISTDP 
  GLOBAL SOFTSTDP 
  GLOBAL EPOTW,EDEPW,IPOTW,IDEPW 
  GLOBAL nextGID 
  GLOBAL EDOPE 
  GLOBAL IDOPE 
  GLOBAL DOPE 
  GLOBAL FORWELIGTR 
  GLOBAL BACKELIGTR 
  GLOBAL EXPELIGTR 
  GLOBAL maxeligtrdur
  GLOBAL reseteligtr 

  
  GLOBAL scaling            
  GLOBAL dynamicdel         
  GLOBAL delspeed           
  GLOBAL scaleinhib         
  GLOBAL activitytau        
  GLOBAL activitybeta       
  GLOBAL activitygamma      
}


PARAMETER {
  tauAM = 10 (ms)
  tauNM = 300 (ms)
  tauGA = 10 (ms)
  tauAM2 = 20 (ms)
  tauNM2 = 300 (ms)
  tauGA2 = 20 (ms)
  invl =  100 (ms)
  WINV =  0
  ahpwt = 0
  tauahp= 10 (ms)
  tauRR = 6 (ms)
  refrac = 5 (ms)
  AHP2REF = 0.0 
  Vbrefrac = 20 (ms)
  RRWght = 0.75
  wwwid = 10
  wwht = 10
  VTH = -45      
  VTHC = -45
  VTHR = -45
  incRR = 0
  Vblock = -20   
  vdt = 0.1      
  mg = 1         
  sop=0
  nbur=1
  tbur=2
  RMP=-65
  EAM = 65
  ENM = 90
  EGA = -15
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
  mg0 = 3.57
  maxnmc = 1.0
  lfpscale = 1.0
  maxplastt = 10.0
  plaststartT = -1 
  plastendT = -1   
  resetplast = 1   
  wsetting = 0 
  ISTDP = 0 
  ESTDP = 1 
  SOFTSTDP = 1 
  EPOTW = 1 
  EDEPW = 1 
  IPOTW = 1
  IDEPW = 1
  nextGID = 0
  DOPE = 0 
  EDOPE = 0 
  IDOPE = 0 
  FORWELIGTR = 1 
  BACKELIGTR = 0 
  EXPELIGTR = 1 
  maxeligtrdur = 100.0 
  reseteligtr = 0 

  
  scaling = 0                          
  dynamicdel = 0                       
  delspeed = 0.0                       
  scaleinhib = 0                       
  activitytau = 100.0e3                
  activitybeta = 4.0e-8                
  activitygamma = 1.0e-10              
}

ASSIGNED {
  Vm VAM VNM VGA AHP VAM2 VNM2 VGA2
  t0 tg twg refractory nxt xloc yloc zloc trrs
  WEX EXSY RES ESIN Psk cbur invlt oinvl tmax spck savclock slowset FLAG
  installed
  pathbeg pathend pathtytarg pathlen
}



CONSTRUCTOR {
  VERBATIM 
  { int lid,lty,lin,lco,lgid,i; unsigned int sz;
    if (ifarg(1)) { lid=(int) *getarg(1); } else { lid= UINT_MAX; } 
    if (ifarg(2)) { lty=(int) *getarg(2); } else { lty= -1; } 
    if (ifarg(3)) { lin=(int) *getarg(3); } else { lin= -1; } 
    if (ifarg(4)) { lco=(int) *getarg(4); } else { lco= -1; } 
    _p_sop = (void*)ecalloc(1, sizeof(id0)); 
    ip = IDP;
    ip->id=lid; ip->type=lty; ip->inhib=lin; ip->col=lco; 
    ip->pg=0x0; ip->dvi=0x0; ip->del=0x0; ip->sprob=0x0; 
    ip->syns=0x0; ip->wgain=0x0; ip->peconv=ip->piconv=0x0; ip->syw1=ip->syw2=0x0;
    ip->pplasttau=0x0; ip->pplastinc=0x0; ip->pplastmaxw=0x0; ip->pdope=0x0;
    ip->dead = ip->invl0 = ip->record = ip->jttr = ip->input = 0; 
    ip->dvt = ip->vbr = ip->wrec = ip->jcn = ip->out = 0;
    for (i=0;i<WRNUM;i++) {ip->wreci[i]=-1; ip->wscale[i]=-1.0;}
    ip->rve=-1;
    pathbeg=-1;
    slowset=0;

    
    
    ip->activity = 0; 
    ip->max_err = 0; 
    ip->max_scale = 100; 
    ip->lastupdate = 0; 
    ip->scalefactor = 1.0; 
    ip->goal_activity = -1; 
    ip->activity_integral_err = 0.0; 
    

    ip->gid = nextGID; nextGID += 1.0;
    process=(int)getpid();
    CNAME[SU]="SU"; CNAME[DP]="DP"; CNAME[IN]="IN";
    if (installed==2.0 && ip->pg) { 
      sz=ivoc_list_count(ip->pg->ce);
      if(verbose) printf("\t**** WARNING new INTF6 created: may want to rerun jitcondiv ****\n");
    } else installed=1.0; 
    cbsv=0x0;
  }
  ENDVERBATIM
}

PROCEDURE resetscaling () {
  VERBATIM
  ip = IDP;
  
  ip->activity = 0; 
  ip->max_err = 0; 
  ip->max_scale = 100; 
  ip->lastupdate = 0; 
  ip->scalefactor = 1.0; 
  ip->goal_activity = -1; 
  ip->activity_integral_err = 0.0; 
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
  trrs = 0
  tmax=0
  pathend=-1
  pathlen=0
  VERBATIM
  { int i,ix;
  ip=IDP;
  _lid=(double)ip->id;
  ip->spkcnt=0;
  ip->blkcnt=0;
  ip->errflag=0;
  ip->pg->lastspk[ip->id]=-1;
  for (i=0;i<CTYN;i++){ix=cty[i]; blockcnt[ix]=spikes[ix]=AMo[ix]=NMo[ix]=GAo[ix]=AMo2[ix]=NMo2[ix]=GAo2[ix]=0;}
  if(seadsetting==3 && resetplast && ip->wgain) for(i=0;i<ip->dvt;i++) ip->wgain[i]=1.0; 
  if(seadsetting==3 && ip->pdope) for(i=0;i<ip->dvt;i++) ip->pdope[i] = -1e9; 
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

  
  
  
  
  
VERBATIM
  activityoneovertau = 1.0 / activitytau; 

ENDVERBATIM
  
}

PROCEDURE reset () {
  Vm = RMP
  VAM = 0
  VNM = 0
  VGA = 0
  AHP=0
  VAM2 = 0
  VNM2 = 0
  VGA2 = 0
  invlt = -1
  t0 = t
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
    sead=hashseed2(2,&x);
  }
  return sead;
}
ENDVERBATIM


FUNCTION DVIDSeed(){
  VERBATIM
  return (double)GetDVIDSeedVal(IDP->id);
  ENDVERBATIM
}


NET_RECEIVE (wAM,wNM,wGA,wGB,wAM2,wNM2,wGA2,wflg) { LOCAL tmp,jcn,id
  INITIAL { wAM=wAM wNM=wNM wGA=wGA wGB=wGB wAM2=wAM2 wNM2=wNM2 wGA2=wGA2 wflg=0}
  
VERBATIM
  id0 *ppre; int prty,poty,prin,prid,poid,ii,sy,nsyn,distal; double STDf,wgain,syw1,syw2; 

ENDVERBATIM
  tmax=t
  VERBATIM
  if (stopoq && !qsz()) stoprun=1;
  ip=IDP; pg=ip->pg; ppre = 0x0; poid=ip->id;
  if (ip->dead) return; 
  _ljcn=ip->jcn; _lid=ip->id;
  tpnt = _pnt; 
  if (PATHMEASURE) { 
    if (_lflag==2 || _lflag<0) { 
      double idty; int i;
      if (_lflag==2) ip->flag=-1; 
      idty=(double)(FOFFSET+ip->id)+1e-2*(double)ip->type+1e-3*(double)ip->inhib+1e-4;
      for (i=0;i<ip->dvt && !stoprun;i++) if (ip->sprob[i]) {
        (*pnt_receive[ip->dvi[i]->_prop->_type])(ip->dvi[i], wts, idty); 
        _p=_pnt->_prop->param; _ppvar=_pnt->_prop->dparam; ip=IDP; 
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

  
  if (dynamicdel) {
    dynamicdelete(t); 
  }
  
  decay_activity_sensor(t); 

  if (scaling) {
    if (ip->goal_activity < 0) {
      
      
      
      
      ip->goal_activity = ip->activity; 
      
    }

    if (!ip->inhib || scaleinhib) {
      
      update_scale_factor(t); 
    }  
  }
  ip->lastupdate = t; 
  


  if (_lflag==OK) { FLAG=OK; flag(); return; } 
  if (_lflag<0) { callback(_lflag); return; }
  pg->eventtot+=1;

  
  ENDVERBATIM
VERBATIM
  if (ip->dbx>2) 
ENDVERBATIM
{ 
    pid() 
    printf("DB0
    if (flag==0) { printf(" (%g %g %g %g %g %g %g)",wAM,wNM,wGA,wAM2,wNM2,wGA2,wflg) }
    printf("\n")
  }


  if (flag>=FOFFSET) { 
    VERBATIM {
      
      prid = (int)(_lflag-FOFFSET); 
      poty=(int)ip->type;
      prty=(int)(1e2*(_lflag-floor(_lflag)));
      prin=(int)(1e3*(_lflag-floor(_lflag)-prty*1e-2)); 
      distal = ((int) (_lflag * 1e5 + 0.5)) % 2;       
      if(distal){ sy=prin?GA2:AM2; } else { sy=prin?GA:AM; }
      
      
      STDf=_args[0]; 
      wgain=_args[1]; 
      syw1=_args[2]; 
      syw2=_args[3]; 
      if(ip->dbx<-1) printf("prid%d,poid%d,wgain=%g\n",prid,poid,wgain); 
      for (ii=0;ii<=6;ii++) _args[ii]=0.; 
      if (seadsetting==3) { 
        ppre = getlp(pg->ce,prid);  
        if(ip->dbx<-1) printf("ppre%p,pre%d->po%d,wg=%g\n",ppre,prid,ip->id,wgain);
        if(ppre->inhib) { 
          if(!ISTDP && !IDOPE) ppre=0x0;
        } else {
          if(!ESTDP && !EDOPE) ppre=0x0;
        }
      }
      if(ppre) { 
        for (ii=sy,nsyn=0;ii<sy+2;ii++) {
          if(ii==AM2 || ii==AM || ii==GA || ii==GA2) { 
            if(wsetting==1.0) { 
              _args[ii] = ii == sy ? syw1 * wgain : syw2 * wgain;
            } else { 
              _args[ii]=wgain*WMAT(prty,poty,ii)*WD0(prty,poty,ii);
            }
            if(ip->dbx<-1) printf("pre%d->po%d,sy=%d,wg=%g,w=%g\n",prid,ip->id,ii,wgain,_args[ii]);
          } else { 
            if(wsetting==1.0) { 
              _args[ii] = ii == sy ? syw1 : syw2;
            } else { 
              _args[ii]=WMAT(prty,poty,ii)*WD0(prty,poty,ii);
            }
          }
          nsyn+=(_args[ii]>0.);
        }
      } else { 
        if(wsetting==1.0) { 
          _args[sy+0] = syw1;
          _args[sy+1] = syw2;
          nsyn = (_args[sy+0]>0.) + (_args[sy+1]>0.);
        } else { 
          for (ii=sy,nsyn=0;ii<sy+2;ii++) nsyn+=((_args[ii]=WMAT(prty,poty,ii)*WD0(prty,poty,ii))>0.);
        }
      }
      if (nsyn==0) return; 

      
      if (scaling) {
        for (ii=sy,nsyn=0;ii<sy+2;ii++) {
          if (!ip->inhib) {
            
            if (ii==AM2 || ii==AM) { 
              
              _args[ii] *= ip->scalefactor;
            }
            if (ii==GA || ii==GA2) {
              
              _args[ii] *= 1 / ip->scalefactor;
            }
          } else {
            
            
            if (ii==AM2 || ii==AM) { 
              
              _args[ii] *= 1 / ip->scalefactor;
            }
            if (ii==GA || ii==GA2) {
              
              _args[ii] *= ip->scalefactor;
            }

          }
        }
      }
      

      if (seadsetting==3) { 
      } else if (seadsetting!=2) { 
        if (seadsetting==1) {
          sead=(unsigned int)(floor(_lflag)*ip->id*seedstep); 
        } else { 
          hsh[0]=floor(_lflag); hsh[1]=(double)ip->id; hsh[2]=seedstep;
          sead=hashseed2(3,&hsh); 
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
      printf(" (%g %g %g %g %g %g %g)",wAM,wNM,wGA,wAM2,wNM2,wGA2,wflg)
      printf("\n")
    }

  } else if (flag==4) { 
    cbur=cbur-1  
    if (cbur>0) { 
      net_send(tbur,4) 
    } else { 
      refractory = 1      
      net_send(refrac-AHP*AHP2REF, 3) 
    }
    tmp=t
VERBATIM
    if (ip->jttr) 
ENDVERBATIM
{ tmp= t+jttr()/10 } 
    if (jcn) { 
      jitcon(tmp)
VERBATIM
      if(ip->out) 
ENDVERBATIM
{ net_event(tmp) } 
    } else { net_event(tmp) }
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
    
    
    
    
  } else if (flag==0 && wflg==1) {
VERBATIM
    ip->input=1; 

ENDVERBATIM
    wflg=2 
    randspk() 
    net_send(nxt,2)
VERBATIM
    return; 

ENDVERBATIM
  } else if (flag==0 && wflg==2) { 
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
  if (VAM2>hoc_epsilon) {VAM2 = VAM2*EXP(-(t - t0)/tauAM2) } else { VAM2=0 } 
  if (VNM2>hoc_epsilon) {VNM2 = VNM2*EXP(-(t - t0)/tauNM2) } else { VNM2=0 } 
  if (VGA2< -hoc_epsilon){VGA2 = VGA2*EXP(-(t - t0)/tauGA2) } else { VGA2=0 } 

  if(refractory==0){
    if(VTHC>VTH) { 
      VTHC = VTH + (VTHR-VTH)*EXP(-(t-trrs)/tauRR) 
    } else if(RRWght<0 && VTHC<VTH) { 
      VTHC = VTH - (VTHR-VTH)*EXP(-(t-trrs)/tauRR) 
    }
  }
  if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } 
  t0 = t 
  Vm = VAM+VNM+VGA+AHP+VAM2+VNM2+VGA2 
  if (Vm> -RMP) {Vm= -RMP}
  if (Vm<  RMP) {Vm= RMP} 

  if (flag==0 || flag>=FOFFSET) { 

    
    if (wAM>0) {
      if (STDAM==0) { VAM = VAM + wAM*(1-Vm/EAM)
      } else        { VAM = VAM + (1-STDAM*STDf)*wAM*(1-Vm/EAM) }
      if (VAM>EAM) { 
VERBATIM
        AMo[ip->type]++; 

ENDVERBATIM
      } else if (VAM<0) { VAM=0 }
    }
    if (wAM2>0) { 
      if (STDAM==0) { VAM2 = VAM2 + wAM2*(1-Vm/EAM)
      } else        { VAM2 = VAM2 + (1-STDAM*STDf)*wAM2*(1-Vm/EAM) }
      if (VAM2>EAM) { 
VERBATIM
        AMo2[ip->type]++; 

ENDVERBATIM
      } else if (VAM2<0) { VAM2=0 }
    }
    
    if (wNM>0 && VNM<ENM) { 
      if (STDNM==0) { VNM = VNM + wNM*rates(RMP+Vm)*(1-Vm/ENM) 
      } else        { VNM = VNM + (1-STDNM*STDf)*wNM*rates(RMP+Vm)*(1-Vm/ENM) }
      if (VNM>ENM) { 
VERBATIM
        NMo[ip->type]++; 

ENDVERBATIM
      } else if (VNM<0) { VNM=0 }
    }
    if (wNM2>0 && VNM2<ENM) { 
      if (STDNM==0) { VNM2 = VNM2 + wNM2*rates(RMP+Vm)*(1-Vm/ENM)
      } else        { VNM2 = VNM2 + (1-STDNM*STDf)*wNM2*rates(RMP+Vm)*(1-Vm/ENM) }
      if (VNM2>ENM) { 
VERBATIM
        NMo2[ip->type]++; 

ENDVERBATIM
      } else if (VNM2<0) { VNM2=0 }
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
          printf("\nAA:%d:%d\n\n",GAo[ip->type],ip->type); 

ENDVERBATIM
        }
      } else if (VGA>0) { VGA=0 } 
    }
    if (wGA2>0 && VGA2>EGA) { 
      if (STDGA==0) {  VGA2 = VGA2 - wGA2*(1-Vm/EGA)
      } else {         VGA2 = VGA2 - (1-STDGA*STDf)*wGA2*(1-Vm/EGA) }
      if (VGA2<EGA) { 
VERBATIM
        GAo2[ip->type]++; 

ENDVERBATIM
VERBATIM
        if (ip->dbx>2) 
ENDVERBATIM
{ 
          pid() printf("DB0A
          if (flag==0) { printf(" (%g %g %g %g %g %g)",wGA2,EGA,VGA2,Vm,AHP,STDf) }  
VERBATIM
          printf("\nAA:%d:%d\n\n",GAo2[ip->type],ip->type); 

ENDVERBATIM
        }
      } else if (VGA2>0) { VGA2=0 } 
    }

VERBATIM
    if (ip->invl0) 
ENDVERBATIM
{ 
      Vm = RMP+VAM+VNM+VGA+AHP+VAM2+VNM2+VGA2
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
      if (jcn) { 
        jitcon(t)
VERBATIM
        if(ip->out) 
ENDVERBATIM
{ net_event(t) } 
      } else { net_event(t) } 
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
  } else if (flag==2) { 
VERBATIM
    if (ip->dbx>1) 
ENDVERBATIM
{pid() printf("DBBa
    if (WEX>1e8) { 
      if (jcn) { 
        jitcon(t)
VERBATIM
        if(ip->out) 
ENDVERBATIM
{ net_event(t) } 
      } else { net_event(t) } 
VERBATIM
      spikes[ip->type]++; 

ENDVERBATIM
      spck=spck+1
VERBATIM
      if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBB
VERBATIM
      if (ip->record) 
ENDVERBATIM
{ recspk(t) } 
VERBATIM
      if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
    } else if (WEX>0) { 
      if(EXSY==AM) {
        tmp = WEX*(1-Vm/EAM)
        VAM = VAM + tmp
      } else if(EXSY==AM2) {
        tmp = WEX*(1-Vm/EAM)
        VAM2 = VAM2 + tmp
      } else if(EXSY==NM) {
        tmp = rates(RMP+Vm)*WEX*(1-Vm/ENM)
        VNM = VNM + tmp
      } else if(EXSY==NM2) {
        tmp = rates(RMP+Vm)*WEX*(1-Vm/ENM)
        VNM2 = VNM2 + tmp
      }
    } else if (WEX<0 && WEX!=-1e9) { 
      if(EXSY==GA) {
        tmp = WEX*(1-Vm/EGA)
        VGA = VGA + tmp
      } else { 
        tmp = WEX*(1-Vm/EGA)
        VGA2 = VGA2 + tmp
      }
    }
    if (WEX!=-1e9) { 
      randspk() 
VERBATIM
      if (ip->input) 
ENDVERBATIM
{ net_send(nxt,2) } 
    }
  } else if (flag==3) { 
    refractory = 0 
    trrs = t 
VERBATIM
    return; 

ENDVERBATIM
  }

  Vm = VAM+VNM+VGA+RMP+AHP+VAM2+VNM2+VGA2 
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
VERBATIM
    

ENDVERBATIM
VERBATIM
    

ENDVERBATIM
VERBATIM
    raise_activity_sensor(t); 

ENDVERBATIM
VERBATIM
    ip->pg->lastspk[ip->id]=_ltmp; 

ENDVERBATIM
VERBATIM
    

ENDVERBATIM
    if (jcn) { 
      jitcon(tmp)
VERBATIM
      if(ip->out) 
ENDVERBATIM
{ net_event(tmp) } 
    } else { net_event(tmp) } 
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
    if(incRR) { 
      VTHC=VTHC+RRWght*(Vblock-VTH)
      if(VTHC > Vblock) {VTHC=Vblock} else if(VTHC < RMP) {VTHC=RMP}
    } else { 
      VTHC=VTH+RRWght*(Vblock-VTH)
    }
    VTHR=VTHC 
    refractory = 1 

    if(seadsetting==3) { 
      if(plaststartT<0 || plastendT<0 || (t>=plaststartT && t<=plastendT)) { 
VERBATIM
        if(ip->dbx<-1) printf("%d@%g applying plasticity\n",ip->id,ip->pg->lastspk[ip->id]); 

ENDVERBATIM
VERBATIM
        if(ESTDP) applyEXSTDP(ip,ip->pg->lastspk[ip->id]); 

ENDVERBATIM
VERBATIM
        if(ISTDP) applyIXSTDP(ip,ip->pg->lastspk[ip->id]); 

ENDVERBATIM
VERBATIM
        if(EDOPE) applyEDOPE(ip,ip->pg->lastspk[ip->id]); 

ENDVERBATIM
VERBATIM
        if(IDOPE) applyIDOPE(ip,ip->pg->lastspk[ip->id]); 

ENDVERBATIM
      }
    }

    if (nbur>1) { 
      cbur=nbur-1 net_send(tbur,4) 
VERBATIM
      return; 

ENDVERBATIM
    } 
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








PROCEDURE jitcon (tm) {
  VERBATIM {
  double mindel, randel, idty, *x; int prty, poty, i, j, k, dv; 
  Point_process *pnt; void* voi;
  
  
  ip=IDP; pg=ip->pg;
  if(verbose>1) printf("col %d , ip %p, pg %p\n",ip->col,ip,pg);
  if (!pg) {printf("No network defined -- must run jitcondiv()\n"); hxe();}
  ip->spkcnt++; 
  if (pg->jrj<pg->jrmax) {  
    pg->jrid[pg->jrj]=(double)ip->id; pg->jrtv[pg->jrj]=_ltm;
    pg->jrj++;
  } else if (wf2 && pg->jrmax) spkoutf2(); 
  pg->jri++;  
  if (jrtm>0) {
    if (t>jrtime) {
      jrtime+=jrtm;
      spkstats2(1.);
    }
  } else if (jrsvd>0 && pg->jri>jrsvn) { 
    jrsvn+=jrsvd; printf("t=%.02f %ld ",t,ip->pg->jri);
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
static int spkoutf2 () {
    fprintf(wf1,"
    fwrite(pg->jrtv,sizeof(double),pg->jrj,wf2); 
    fwrite(pg->jrid,sizeof(double),pg->jrj,wf2); 
    fflush(wf1); fflush(wf2);
    pg->jrj=0;
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
  ip=IDP; pg=ip->pg;
  flag=(int)(_lflag+1e-6);
  clk=clock()-savclock; savclock=clock();
  if (cbsv) hoc_call_func(cbsv,0);
  if (tf) fprintf(tf,"t=%.02f;%ld(%g) ",t,pg->jri,clk/1e6); else {
    printf("t=%.02f;%ld(%g) ",t,pg->jri,clk/1e6); }
  for (i=0;i<CTYN;i++) {
    ix=cty[i];
    pg->spktot+=spikes[ix];
    if (tf) {
      fprintf(tf,"%s:%d/%d:%d;%d;%d;%d;%d;%d ",CNAME[i],spikes[ix],\
              blockcnt[ix],AMo[ix],NMo[ix],GAo[ix],AMo2[ix],NMo2[ix],GAo2[ix]);
    } else {
      printf("%s:%d/%d:%d;%d;%d;%d;%d;%d ",CNAME[i],spikes[ix],blockcnt[ix],\
             AMo[ix],NMo[ix],GAo[ix],AMo2[ix],NMo2[ix],GAo2[ix]);
    }
    spck=0;
    blockcnt[ix]=spikes[ix]=0;
    AMo[ix]=NMo[ix]=GAo[ix]=AMo2[ix]=NMo2[ix]=GAo2[ix]=0;
  }
  if (tf && flag==2) {  fprintf(tf,"\nt=%g tot_spks: %ld; tot_events: %ld\n",t,pg->spktot,pg->eventtot); 
  } else if (flag==2) {  printf("\ntotal spikes: %ld; total events: %ld\n",pg->spktot,pg->eventtot); 
  } else if (tf) fprintf(tf,"\n"); else printf("\n");
}
ENDVERBATIM
}

PROCEDURE oobpr () {
VERBATIM {
  int i,ix;
  for (i=0;i<CTYN;i++){ 
    ix=cty[i];
    printf("%d:%d/%d:%d;%d;%d;%d;%d;%d ",ix,spikes[ix],blockcnt[ix],\
           AMo[ix],NMo[ix],GAo[ix],AMo2[ix],NMo2[ix],GAo2[ix]);
  }
  printf("\n");
}
ENDVERBATIM
}

PROCEDURE callback (fl) {
  VERBATIM {
  int i; double idty, idtflg, del0, ddel; id0 *jp; Point_process *upnt; 
  i=(unsigned int)((-_lfl)-1); 
  jp=IDP; upnt=tpnt; del0=jp->del[i]; ddel=0.;
  idty=(double)(FOFFSET+jp->id)+1e-2*(double)jp->type+1e-3*(double)jp->inhib+1e-4;
  while (ddel<=DELMIN) { 
    if (Vblock<VTHC) { 
      wts[0]=0; 
    } else if(STDAM || STDNM || STDGA) { 
                                                         
                                                         
                                                         
      wts[0]=(VTHC-VTH)/(Vblock-VTH); 
    }
    wts[1]=0.0; 
    if(seadsetting==3) { 
      if(jp->inhib) {
        if(ISTDP || IDOPE) wts[1]=jp->wgain[i];
      } else {
        if(ESTDP || EDOPE) wts[1]=jp->wgain[i];
      }
    }
    if(wsetting==1.0 && jp->syw1 && jp->syw2) {wts[2]=jp->syw1[i]; wts[3]=jp->syw2[i]; } 
    idtflg = idty + (1e-5 * jp->syns[i]);
    
    if (jp->sprob[i]) (*pnt_receive[jp->dvi[i]->_prop->_type])(jp->dvi[i], wts, idtflg); 
    _p=upnt->_prop->param; _ppvar=upnt->_prop->dparam; 
    i++;
    if (i>=jp->dvt) return 0; 
    ddel=jp->del[i]-del0;   
  }
  
  while (i<jp->dvt && (!jp->sprob[i] || (*(id0**)&(jp->dvi[i]->_prop->dparam[2]))->dead)) i++;
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
        if (!(lb=ivoc_list_item(pg->ce,(unsigned int)floor(dscr[j]+pg->ix[poty])))) {
          printf("INTF6:callback %g exceeds %d for list ce\n",floor(dscr[j]+pg->ix[poty]),pg->cesz); 
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
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt; ip->syns=(char*)calloc(dvt,sizeof(char));
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); 
  for (i=0;i<dvt;i++) ip->sprob[i]=1; 
  free(da); free(db); 
  }
ENDVERBATIM
}


PROCEDURE patha2b () {
  VERBATIM
  int i; double idty, *x; static Point_process *_pnt; static id0 *ip0;
  ip=IDP; pg=ip->pg;
  pathbeg=*getarg(1); pathidtarg=*getarg(2);
  pathtytarg=-1;  PATHMEASURE=1; pathlen=stopoq=0;
  for (i=0;i<pg->cesz;i++) { lop(pg->ce,i); 
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
  ip=IDP; pg=ip->pg;
  x=0x0;
  s=hoc_lookup("finitialize");
  if (ifarg(2)) {
    na=vector_arg_px(1,&a);
    nb=vector_arg_px(2,&b);
    if (ifarg(3)) x=vector_newsize(vector_arg(3),na*nb);
  } else {
    na=nb=pg->cesz;  
    if (ifarg(1)) x=vector_newsize(vector_arg(1),na*nb);
  }
  
  pfl = (char **)malloc(pg->cesz * (unsigned)sizeof(char *));
  for (i=0;i<pg->cesz;i++) { lop(pg->ce,i); scr[i]=qp->inhib; pfl[i]=&qp->flag; }
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
      
      for (i=0;i<pg->cesz;i++) *pfl[i]=0;
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
    int i,j,k,iarg,av1,a2,a3,a4,a6,a7,dvt,getactive=0,idx=0,*pact,prty,poty,sy,ii; 
    double *dbs, *x,*x1,*x2,*x3,*x4,*x5,*x6,*x7,idty,y[2],flag;
    void* voi, *voi2,*voi3; Point_process **das;
    ip=IDP; pg=ip->pg;
    getactive=a2=a3=a4=0;
    if (ip->dead) return 0.0;
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
    if (flag==2) { x1=vector_newsize(voi,CTYPi); for (i=0;i<CTYPi;i++) x1[i]=0;
    } else x1=vector_newsize(voi,dvt);
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x2=vector_newsize(voi,dvt);  a2=1; }
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x3=vector_newsize(voi,dvt); a3=1;}
    if (ifarg(iarg)) { 
      voi=vector_arg(iarg++); x4=vector_newsize(voi,dvt); a4=1;
      voi=vector_arg(iarg++); x5=vector_newsize(voi,dvt);
    }
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x6=vector_newsize(voi,dvt); a6=1;} else a6=0;
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x7=vector_newsize(voi,dvt); a7=1;} else a7=0;
    idty=(double)(FOFFSET+ip->id)+1e-2*(double)ip->type+1e-3*(double)ip->inhib+1e-4;
    prty=ip->type; sy=ip->inhib?GA:AM;
    for (i=0,j=0;i<dvt;i++) {
      qp=*((id0**) &((das[i]->_prop->dparam)[2])); 
      if (getactive && (qp->dead || ip->sprob[i]==0)) continue;
      if (flag==1) { x1[j]=(double)qp->type; 
      } else if (flag==2) { x1[qp->type]++; 
      } else if (flag==3) { x1[j]=(double)qp->col; 
      } else x1[j]=(double)qp->id;
      if (a2) x2[j]=dbs[i];
      if (a3) x3[j]=(double)ip->sprob[i];
      if (a4) {
        if(ip->inhib){sy=ip->syns[i]?GA2:GA;} else {sy=ip->syns[i]?AM2:AM;} 
        poty = qp->type;
        if(wsetting==1) { 
          y[0]=ip->syw1[i]; y[1]=ip->syw2[i];
        } else {
          if (seadsetting==2 || seadsetting==3) { 
            for(ii=0;ii<2;ii++) y[ii]=WMAT(prty,poty,sy+ii)*WD0(prty,poty,sy+ii);
          } else {
            if (seadsetting==1) { 
              sead=(unsigned int)(FOFFSET+ip->id)*qp->id*seedstep; 
            } else { 
              hsh[0]=(double)(FOFFSET+ip->id); hsh[1]=(double)(qp->id); hsh[2]=seedstep;
              sead=hashseed2(3,&hsh); 
            }
            mcell_ran4(&sead, y, 2, 1.);
            for(ii=0;ii<2;ii++) {
              y[ii]=2*WVAR*(y[ii]+0.5/WVAR-0.5)*WMAT(prty,poty,sy+ii)*WD0(prty,poty,sy+ii); }
          }
        }
        x4[j]=y[0]; x5[j]=y[1];
      }
      if (a6) x6[j] = ip->syns[i];  
      if (a7 && ip->wgain)x7[j]=ip->wgain[i];
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
  void* voi; Point_process **das; id0 *pp;
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
  for (i=0,k=0; i<pg->cesz; i++) {
    lop(pg->ce,i);
    if (getactive && qp->dead) continue;
    dvt=qp->dvt; das=qp->dvi;
    for (j=0;j<dvt;j++) {
      if (getactive && qp->sprob[j]==0) continue;
      if (ip==*((id0**) &((das[j]->_prop->dparam)[2]))) {
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
  ip=IDP; pg=ip->pg;
  int iListSz=ivoc_list_count(pList),iCell,iStartID=ifarg(2)?*getarg(2):0,\
    iEndID=ifarg(3)?*getarg(3):pg->cesz-1;
  int skipinhib = ifarg(4)?*getarg(4):0, i,j,nv,*pused=(int*)calloc(pg->cesz,sizeof(int)),iSyns=0;
  double **vvo = (double**)malloc(sizeof(double*)*iListSz),\
    *psyns=(double*)calloc(pg->cesz,sizeof(double));
  id0* rp;
  for(iCell=iStartID;iCell<=iEndID;iCell++){
    if(verbose && iCell%1000==0) printf("%d ",iCell);
    lop(pg->ce,iCell);
    if(!qp->dvt || (skipinhib && qp->inhib)){
      list_vector_resize(pList,iCell,0);
      continue;
    }
    iSyns=0;
    for(j=0;j<qp->dvt;j++){      
      rp=*((id0**) &((qp->dvi[j]->_prop->dparam)[2])); 
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
  ip=IDP; pg=ip->pg;
  printf("reading: ");
  for(iCell=0;iCell<pg->cesz;iCell++){
    if(iCell%1000==0)printf("%d ",iCell);
    lop(pg->ce,iCell);
    fread(&qp->id,sizeof(unsigned int),1,fp); 
    fread(&qp->type,sizeof(unsigned char),1,fp); 
    fread(&qp->col,sizeof(unsigned int),1,fp); 
    fread(&qp->dead,sizeof(unsigned char),1,fp); 
    fread(&qp->dvt,sizeof(unsigned int),1,fp); 
    
    if(qp->del){ free(qp->del); free(qp->dvi); free(qp->sprob);
      qp->dvt=0; qp->dvi=(Point_process**)0x0; qp->del=(double*)0x0; qp->sprob=(char *)0x0; }
    
    if(!qp->dvt) continue;
    qp->dvi = (Point_process**)malloc(sizeof(Point_process*)*qp->dvt);  
    for(i=0;i<qp->dvt;i++){
      fread(&iOutID,sizeof(unsigned int),1,fp); 
      if (!(lb=ivoc_list_item(pg->ce,iOutID))) {
        printf("INTF6:callback %d exceeds %d for list ce\n",iOutID,pg->cesz); hxe(); }
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
  ip=IDP; pg=ip->pg;
  printf("writing: ");
  for(iCell=0;iCell<pg->cesz;iCell++){
    if(iCell%1000==0)printf("%d ",iCell);
    lop(pg->ce,iCell);
    fwrite(&qp->id,sizeof(unsigned int),1,fp); 
    fwrite(&qp->type,sizeof(unsigned char),1,fp); 
    fwrite(&qp->col,sizeof(unsigned int),1,fp); 
    fwrite(&qp->dead,sizeof(unsigned char),1,fp); 
    fwrite(&qp->dvt,sizeof(unsigned int),1,fp); 
    if(!qp->dvt)continue; 
    for(i=0;i<qp->dvt;i++){
      pnnt=qp->dvi[i];
      fwrite(&(*(id0**)&(pnnt->_prop->dparam[2]))->id,sizeof(unsigned int),1,fp); 
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
  ip=IDP; pg=ip->pg;
  pListWires = AllocListVec(*hoc_objgetarg(1));
  idvfl=flag=0; iStartID=0; iEndID=pg->cesz-1;
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
    lop(pg->ce,iCell);
    if (qp->dead) continue;
    y=pListWires->pv[i]; dvt=pListWires->plen[i];
    if(!dvt) continue; 
    d=pListDels->pv[i];  dn=pListDels->plen[i];
    if (dn!=dvt) {printf("setdvir() ERR vec sizes for wire,delay list entries not equal %d: %d %d\n",i,dvt,dn); hxe();}
    setdvi2(y,d,0x0,dvt,flag,0x0,0x0);
  }
  FreeListVec(&pListWires);
  FreeListVec(&pListDels);
  return 1.0;
  ENDVERBATIM
}

PROCEDURE clrdvi () {
  VERBATIM
  int i;
  ip=IDP; pg=ip->pg;
  for (i=0;i<pg->cesz;i++) { 
    lop(pg->ce,i);
    if (qp->dvt!=0x0) {
      free(qp->dvi); free(qp->del); free(qp->sprob);
      qp->dvt=0; qp->dvi=(Point_process**)0x0; qp->del=(double*)0x0; qp->sprob=(char *)0x0;
      if(wsetting==1) freesywv(qp);
    }
  }
  ENDVERBATIM
}


PROCEDURE setdviv () {
  VERBATIM
  int i,j,k,l,nprv,dvt,*scr; double *prv,*pov,*dlv,x,*ds,*w1,*w2; char* s;
  ip=IDP; pg=ip->pg;
  nprv=vector_arg_px(1, &prv);
  i=vector_arg_px(2, &pov);
  j=vector_arg_px(3, &dlv);
  if(ifarg(4)) { s=(char*)calloc((l=vector_arg_px(4,&ds)),sizeof(char)); for(k=0;k<l;k++) s[k]=ds[k]; k=0;
  } else s=0x0;
  if (nprv!=i || i!=j || j!=l) {printf("intf:setdviv ERRA: %d %d %d %d\n",nprv,i,j,l); hxe();}
  if (wsetting==1) {
    i=vector_arg_px(5, &w1);
    j=vector_arg_px(6, &w2);
    if (nprv!=i || i!=j) {printf("intf:setdviv ERRB: %d %d %d\n",nprv,i,j); hxe();}
  }
  
  scr=(int *)ecalloc(pg->cesz, sizeof(int));
  for (i=0;i<pg->cesz;i++) scr[i]=0;
  for (i=0,j=-1;i<nprv;i++) {
    if ((int)prv[i]<j) { printf("intf:setdviv ERRC vecs should be sorted by prid vec\n");hxe(); }
    j=(int)prv[i];
    scr[j]++;
  }
  if (ip->dbx>1) for (i=0;i<pg->cesz;i++) printf("%d ",scr[i]);
  for (i=-1,k=0;k<nprv;k+=dvt) { if(i%1000==0) printf(".");
    if ((int)prv[k]==i) {printf("intf:setdviv ERRD number repeated %g %d %d\n",prv[k],i,k);hxe();}
    i=(int)prv[k]; 
    lop(pg->ce,i); 
    dvt=scr[i]; 
    if (ip->dbx>0) printf("DBA:%d,%d,%d ",i,dvt,k);
    if (qp->dead) continue;
    if (dvt>0) {
      if (wsetting==1) {
        setdvi3(pov+k,dlv+k,s+k,dvt,w1+k,w2+k); 
      } else {
        setdvi2(pov+k,dlv+k,s?s+k:0x0,dvt,1,0x0,0x0);
      }
    }
  }
  if(s) free(s);
  ENDVERBATIM
}

VERBATIM
void setupsywv (id0* p, int sz) {
  p->syw1 = p->syw1!=0x0 ? (double*) realloc((double*) p->syw1, sz*sizeof(double)) :
                           (double*) calloc(sz,sizeof(double));

  p->syw2 = p->syw2!=0x0 ? (double*) realloc((double*) p->syw2, sz*sizeof(double)) :
                           (double*) calloc(sz,sizeof(double));
}








void freesywv (id0* p) {
  if(p->syw1) free(p->syw1); p->syw1=0x0;
  if(p->syw2) free(p->syw2); p->syw2=0x0;
}
ENDVERBATIM


FUNCTION setsywv () {
  VERBATIM
  int sz,n1,n2; double *psyw1,*psyw2; id0* ip;
  ip=IDP; pg=ip->pg; sz=ip->dvt;
  if((n1=vector_arg_px(1, &psyw1))!=sz || (n2=vector_arg_px(2, &psyw2))!=sz) {
    printf("setsywv ERRA: make sure weight vector sizes (%d,%d) same size as div(%d)\n",n1,n2,sz);
    return 0.0;
  }
  setupsywv(ip,sz); 
  memcpy(ip->syw1,psyw1,sizeof(double)*sz); 
  memcpy(ip->syw2,psyw2,sizeof(double)*sz);
  return sz;
  ENDVERBATIM
}


FUNCTION getsywv () {
  VERBATIM
  int sz,n1,n2; double *psyw1,*psyw2; id0* ip;
  ip=IDP; pg=ip->pg; sz=ip->dvt;
  if(!ip->syw1 || !ip->syw2) {
    printf("getsywv ERRA: syw1,syw2 were never initialized with setsywv!\n");
    return 0.0;
  }
  if((n1=vector_arg_px(1, &psyw1))!=sz || (n2=vector_arg_px(2, &psyw2))!=sz) {
    printf("getsywv ERRB: make sure weight vector sizes (%d,%d) same size as div(%d)\n",n1,n2,sz);
    return 0.0;
  }
  memcpy(psyw1,ip->syw1,sizeof(double)*sz); 
  memcpy(psyw2,ip->syw2,sizeof(double)*sz);
  return sz;
  ENDVERBATIM
}

VERBATIM

int* getpeconv (id0* ip,int* psz) {
  Point_process **das; int* pfrom;
  int i,j,k,dvt;
  *psz=ip->dvt>0?ip->dvt:16; pg=ip->pg;
  pfrom=(int*) calloc(psz[0],sizeof(int));
  for (i=0,k=0; i<pg->cesz; i++) {
    lop(pg->ce,i);
    if(qp->inhib) continue; 
    dvt=qp->dvt;
    das=qp->dvi;
    for (j=0;j<dvt;j++) {
      if (ip==*((id0**) &((das[j]->_prop->dparam)[2]))) {
        if (k>=*psz) {
          psz[0]*=2;
          pfrom=(int*) realloc((void*)pfrom,psz[0]*sizeof(int));
        }
        pfrom[k]=qp->id;
        k++;
        break;
      }
    }
  }
  *psz=k;
  return pfrom;
}


int* getpiconv (id0* ip,int* psz) {
  Point_process **das; int* pfrom;
  int i,j,k,dvt;
  *psz=ip->dvt>0?ip->dvt:16; pg=ip->pg;
  pfrom=(int*) calloc(psz[0],sizeof(int));
  for (i=0,k=0; i<pg->cesz; i++) {
    lop(pg->ce,i);
    if(!qp->inhib) continue; 
    dvt=qp->dvt;
    das=qp->dvi;
    for (j=0;j<dvt;j++) {
      if (ip==*((id0**) &((das[j]->_prop->dparam)[2]))) {
        if (k>=*psz) {
          psz[0]*=2;
          pfrom=(int*) realloc((void*)pfrom,psz[0]*sizeof(int));
        }
        pfrom[k]=qp->id;
        k++;
        break;
      }
    }
  }
  *psz=k;
  return pfrom;
}


int myfindidx (id0* ppre,int poid) {
  int i; Point_process** das; id0* ppo;
  das=ppre->dvi;
  for(i=0;i<ppre->dvt;i++) {
    ppo=*((id0**) &((das[i]->_prop->dparam)[2])); 
    if(ppo->id==poid) return i;
  }
  return -1;
}



static void applyEDOPE (id0* pcell,double myspkt) {
  int poid,prid,sz,i,idx; postgrp* pg; double d,inc,tmp,pinc,tau,maxw; id0* ppre, *ppo;
  if(seadsetting!=3.) return; 
  poid=pcell->id; pg=pcell->pg;
  if(pcell->dbx<-1) printf("applyEDOPE: pcell=%p\n",pcell);
  if (FORWELIGTR) {  
    for(i=0;i<pcell->econvsz;i++) {
      prid = pcell->peconv[i];                
      if(pg->lastspk[prid]<0) continue;     
      if( (d = myspkt - pg->lastspk[prid] ) > maxplastt) continue;  
      if(verbose>2) printf("spk%d:%g, spk%d:%g, d=%g\n",prid,pg->lastspk[prid],poid,pg->lastspk[poid],d);
      ppre = getlp(pg->ce,prid);            
      idx = myfindidx(ppre,poid);           
      if(idx<0){printf("**** applyEDOPE ERR: bad idx = %d!!!!!!!!!\n",idx); return;}
      if( ! ( inc = ppre->pplastinc[idx] ) ) continue; 
      ppre->pdope[idx] = t; 
      if(verbose>2) printf("EDOPEA:ppre->inhib=%d,pcell->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d\n",ppre->inhib,pcell->inhib,d,tau,d/tau,prid,poid);
    }
  }
  if (BACKELIGTR) { 
    if(pcell->inhib) return; 
    for(i=0;i<pcell->dvt;i++) { 
      ppo=*((id0**) &((pcell->dvi[i]->_prop->dparam)[2])); 
      poid = ppo->id;
      if(pg->lastspk[poid]<0) continue;
      if( (d = myspkt - pg->lastspk[poid] ) <= maxplastt) {
        if( ! ( inc = pcell->pplastinc[i] ) ) continue;
        pcell->pdope[i] = -t; 
        if(verbose>2) printf("EDOPEB:ppo->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d\n",ppo->inhib,d,tau,d/tau,prid,poid);
      }    
    } 
  }
}






static void applyIDOPE (id0* pcell,double myspkt) {
  int poid,prid,sz,i,idx; postgrp* pg; double d,inc,tmp,pinc,tau,maxw; id0* ppre, *ppo;
  if(seadsetting!=3.) return; 
  poid=pcell->id; pg=pcell->pg;
  if(pcell->dbx<-1) printf("applyplast: pcell=%p\n",pcell);
  if (FORWELIGTR) {  
    for(i=0;i<pcell->iconvsz;i++) {
      prid = pcell->piconv[i];                
      if(pg->lastspk[prid]<0) continue;     
      if( (d = myspkt - pg->lastspk[prid] ) > maxplastt) continue;  
      if(verbose>2) printf("spk%d:%g, spk%d:%g, d=%g\n",prid,pg->lastspk[prid],poid,pg->lastspk[poid],d);
      ppre = getlp(pg->ce,prid);            
      idx = myfindidx(ppre,poid);           
      if(idx<0){printf("**** applyISSTDP ERR: bad idx = %d!!!!!!!!!\n",idx); return;}
      if( ! ( inc = ppre->pplastinc[idx] ) ) continue; 
      ppre->pdope[idx] = t; 
      if(verbose>2) printf("IDOPEA:ppre->inhib=%d,pcell->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d\n",ppre->inhib,pcell->inhib,d,tau,d/tau,prid,poid);
    }
  }
  if(BACKELIGTR) { 
    if(!pcell->inhib) return; 
    for(i=0;i<pcell->dvt;i++) { 
      ppo=*((id0**) &((pcell->dvi[i]->_prop->dparam)[2])); 
      poid = ppo->id;
      if(pg->lastspk[poid]<0) continue;
      if( (d = myspkt - pg->lastspk[poid] ) <= maxplastt) {
        if( ! ( inc = pcell->pplastinc[i] ) ) continue;
        pcell->pdope[i] = -t; 
        if(verbose>2) printf("IDOPEB:ppo->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d\n",ppo->inhib,d,tau,d/tau,prid,poid);
      }    
    }
  }
}



static void applyEXSTDP (id0* pcell,double myspkt) {
  int poid,prid,sz,i,idx; postgrp* pg; double d,inc,tmp,pinc,tau,maxw; id0* ppre, *ppo;
  if(seadsetting!=3.) return; 
  poid=pcell->id; pg=pcell->pg;
  if(pcell->dbx<-1) printf("applyEXSTDP: pcell=%p\n",pcell);
  for(i=0;i<pcell->econvsz;i++) {
    prid = pcell->peconv[i];                
    if(pg->lastspk[prid]<0) continue;     
    if( (d = myspkt - pg->lastspk[prid] ) > maxplastt) continue;  
    if(verbose>2) printf("spk%d:%g, spk%d:%g, d=%g\n",prid,pg->lastspk[prid],poid,pg->lastspk[poid],d);
    ppre = getlp(pg->ce,prid);            
    idx = myfindidx(ppre,poid);           
    if(idx<0){printf("**** applyEXSTDP ERR: bad idx = %d!!!!!!!!!\n",idx); return;}
    if( ! ( inc = ppre->pplastinc[idx] ) ) continue; 
    tau = ppre->pplasttau[idx];
    maxw = ppre->pplastmaxw[idx];
    tmp = ppre->wgain[idx]; 
    if(SOFTSTDP) inc *= (1.0 - tmp / maxw); 
    ppre->wgain[idx] += EPOTW * inc * exp( -d / tau ); 
    if(ppre->wgain[idx]<0.) ppre->wgain[idx]=0.; 
    else if(!SOFTSTDP && ppre->wgain[idx]>maxw) ppre->wgain[idx]=maxw;
    if(verbose>2) printf("PLAST:ppre->inhib=%d,pcell->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d: inc=%g, wgA=%g, wgB=%g\n",ppre->inhib,pcell->inhib,d,tau,d/tau,prid,poid,inc,tmp,ppre->wgain[idx]);
  }
  if(pcell->inhib) return; 
  for(i=0;i<pcell->dvt;i++) { 
    ppo=*((id0**) &((pcell->dvi[i]->_prop->dparam)[2])); 
    poid = ppo->id;
    if(pg->lastspk[poid]<0) continue;
    if( (d = myspkt - pg->lastspk[poid] ) < maxplastt) {
      if( ! ( inc = pcell->pplastinc[i] ) ) continue;
      tau = pcell->pplasttau[i];
      maxw = pcell->pplastmaxw[i];
      tmp = pcell->wgain[i]; 
      if(SOFTSTDP) inc *= (tmp / maxw); 
      pcell->wgain[i] -= EDEPW * inc * exp( -d / tau ); 
      if(pcell->wgain[i]<0.) pcell->wgain[i]=0.; 
      else if(!SOFTSTDP && pcell->wgain[i]>maxw) pcell->wgain[i]=maxw;
      if(verbose>2) printf("DEP:ppo->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d: inc=%g, wgA=%g, wgB=%g\n",ppo->inhib,d,tau,d/tau,prid,poid,inc,tmp,pcell->wgain[i]);
    }    
  }
}



static void applyIXSTDP (id0* pcell,double myspkt) {
  int poid,prid,sz,i,idx; postgrp* pg; double d,inc,tmp,pinc,tau,maxw; id0* ppre, *ppo;
  if(seadsetting!=3.) return; 
  poid=pcell->id; pg=pcell->pg;
  if(pcell->dbx<-1) printf("applyplast: pcell=%p\n",pcell);
  for(i=0;i<pcell->iconvsz;i++) {
    prid = pcell->piconv[i];                
    if(pg->lastspk[prid]<0) continue;     
    if( (d = myspkt - pg->lastspk[prid] ) > maxplastt) continue;  
    if(verbose>2) printf("spk%d:%g, spk%d:%g, d=%g\n",prid,pg->lastspk[prid],poid,pg->lastspk[poid],d);
    ppre = getlp(pg->ce,prid);            
    idx = myfindidx(ppre,poid);           
    if(idx<0){printf("**** applyISSTDP ERR: bad idx = %d!!!!!!!!!\n",idx); return;}
    if( ! ( inc = ppre->pplastinc[idx] ) ) continue; 
    tau = ppre->pplasttau[idx];
    maxw = ppre->pplastmaxw[idx];
    tmp = ppre->wgain[idx]; 
    if(SOFTSTDP) inc *= (tmp / maxw); 
    ppre->wgain[idx] -= IDEPW * inc * exp( -d / tau ); 
    if(ppre->wgain[idx]<0.) ppre->wgain[idx]=0.; 
    else if(!SOFTSTDP && ppre->wgain[idx]>maxw) ppre->wgain[idx]=maxw;
    if(verbose>2) printf("DEP:ppre->inhib=%d,pcell->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d: inc=%g, wgA=%g, wgB=%g\n",ppre->inhib,pcell->inhib,d,tau,d/tau,prid,poid,inc,tmp,ppre->wgain[idx]);
  }
  if(!pcell->inhib) return; 
  for(i=0;i<pcell->dvt;i++) { 
    ppo=*((id0**) &((pcell->dvi[i]->_prop->dparam)[2])); 
    poid = ppo->id;
    if(pg->lastspk[poid]<0) continue;
    if( (d = myspkt - pg->lastspk[poid] ) < maxplastt) {
      if( ! ( inc = pcell->pplastinc[i] ) ) continue;
      tau = pcell->pplasttau[i];
      maxw = pcell->pplastmaxw[i];
      tmp = pcell->wgain[i]; 
      if(SOFTSTDP) inc *= (1.0 - tmp / maxw); 
      pcell->wgain[i] += IPOTW * inc * exp( -d / tau ); 
      if(pcell->wgain[i]<0.) pcell->wgain[i]=0.; 
      else if(!SOFTSTDP && pcell->wgain[i]>maxw) pcell->wgain[i]=maxw;
      if(verbose>2) printf("PLAST:ppo->inhib=%d,d=%g,tau=%g,d/tau=%g, %d->%d: inc=%g, wgA=%g, wgB=%g\n",ppo->inhib,d,tau,d/tau,prid,poid,inc,tmp,pcell->wgain[i]);
    }    
  }
}

ENDVERBATIM


FUNCTION geteconv () {
  VERBATIM
  int i; double *x; void *voi;
  ip=IDP; pg=ip->pg;
  if(!ip->peconv) ip->peconv=getpeconv(ip,&ip->econvsz);
  voi=vector_arg(1);
  x=vector_newsize(voi,ip->econvsz);
  for(i=0;i<ip->econvsz;i++) x[i]=(double)ip->peconv[i];
  return ip->econvsz;
  ENDVERBATIM
}


FUNCTION geticonv () {
  VERBATIM
  int i; double *x; void *voi;
  ip=IDP; pg=ip->pg;
  if(!ip->piconv) ip->piconv=getpiconv(ip,&ip->iconvsz);
  voi=vector_arg(1);
  x=vector_newsize(voi,ip->iconvsz);
  for(i=0;i<ip->iconvsz;i++) x[i]=(double)ip->piconv[i];
  return ip->iconvsz;
  ENDVERBATIM
}


VERBATIM
static int finishdvi2 (struct ID0* p) {
  Point_process **da,**das;
  double *db,*dbs,*w1,*w1s,*w2,*w2s;
  char *syns,*synss;
  int i, dvt;
  db=p->del;
  da=p->dvi; 
  dvt=p->dvt;
  syns=p->syns;
  dbs=(double*)malloc(dvt*sizeof(double)); 
  das=(Point_process**)malloc(dvt*sizeof(Point_process*)); 
  synss=(char*)malloc(dvt*sizeof(char)); 
  if(wsetting==1 && p->syw1 && p->syw2) {
    w1=p->syw1;
    w2=p->syw2;
    w1s=(double*)malloc(dvt*sizeof(double)); 
    w2s=(double*)malloc(dvt*sizeof(double));
    gsort5(db,da,syns,w1,w2,dvt,dbs,das,synss,w1s,w2s); 
    p->syw1=w1s; p->syw2=w2s; 
    free(w1); free(w2); 
  } else gsort3(db,da,syns,dvt,dbs,das,synss);
  p->del=dbs; p->dvi=das; p->syns=synss;
  free(db); free(da); free(syns); 
  p->sprob=(unsigned char*)realloc((void*)p->sprob,(size_t)dvt*sizeof(char));
  for (i=0;i<dvt;i++) p->sprob[i]=1; 
  p->wgain=(double*)realloc((void*)p->wgain,(size_t)dvt*sizeof(double));
  for (i=0;i<dvt;i++) p->wgain[i]=1.0; 
  p->peconv = getpeconv(p,&p->econvsz); 
  p->piconv = getpiconv(p,&p->iconvsz); 
  if(seadsetting==3) {
    p->pplasttau = (double*)realloc((void*)p->pplasttau,(size_t)dvt*sizeof(double));
    p->pplastinc = (double*)realloc((void*)p->pplastinc,(size_t)dvt*sizeof(double));
    p->pplastmaxw = (double*)realloc((void*)p->pplastmaxw,(size_t)dvt*sizeof(double));
    if(DOPE) p->pdope = (double*)realloc((void*)p->pdope,(size_t)dvt*sizeof(double));
  }
}
ENDVERBATIM


PROCEDURE finishdvir () {
  VERBATIM
  int iCell;
  ip=IDP; pg=ip->pg;
  for(iCell=0;iCell<pg->cesz;iCell++){
    lop(pg->ce,iCell);
    finishdvi2(qp);
  }
  ENDVERBATIM
}


PROCEDURE finishdvi () {
VERBATIM
  finishdvi2(IDP);
ENDVERBATIM
}



FUNCTION setplast () {
  VERBATIM
  double *wgain,*pplasttau,*pplastinc,*pplastmaxw;
  if(seadsetting!=3) {printf("setplast ERR0: seadsetting must be 3, plast mode off!\n"); return 0;}
  ip=IDP; pg=ip->pg;
  if(vector_arg_px(1,&wgain) != ip->dvt ||
     vector_arg_px(2,&pplasttau) != ip->dvt ||
     vector_arg_px(3,&pplastinc) != ip->dvt ||
     vector_arg_px(4,&pplastmaxw) != ip->dvt) {printf("setplast ERR1: input vectors must have size %d!\n",ip->dvt); return 0;}
  memcpy(ip->wgain,wgain,sizeof(double)*ip->dvt);
  memcpy(ip->pplasttau,pplasttau,sizeof(double)*ip->dvt);
  memcpy(ip->pplastinc,pplastinc,sizeof(double)*ip->dvt);
  memcpy(ip->pplastmaxw,pplastmaxw,sizeof(double)*ip->dvt);
  return 0.0;
  ENDVERBATIM
}



FUNCTION getplast () {
  VERBATIM
  double *wgain,*pplasttau,*pplastinc,*pplastmaxw;
  if(seadsetting!=3) {printf("getplast ERR0: seadsetting must be 3, plast mode off!\n"); return 0;}
  ip=IDP; pg=ip->pg;
  ip=IDP; pg=ip->pg;
  if(vector_arg_px(1,&wgain) != ip->dvt ||
     vector_arg_px(2,&pplasttau) != ip->dvt ||
     vector_arg_px(3,&pplastinc) != ip->dvt ||
     vector_arg_px(4,&pplastmaxw) != ip->dvt) {printf("getplast ERR1: output vectors must have size %d!\n",ip->dvt); return 0;}
  memcpy(wgain,ip->wgain,sizeof(double)*ip->dvt);
  memcpy(pplasttau,ip->pplasttau,sizeof(double)*ip->dvt);
  memcpy(pplastinc,ip->pplastinc,sizeof(double)*ip->dvt);
  memcpy(pplastmaxw,ip->pplastmaxw,sizeof(double)*ip->dvt);
  return 1.0;
  ENDVERBATIM
}



PROCEDURE setdvi () {
VERBATIM {
  int i,j,k,dvt,flag; double *d, *y, *ds, *w1, *w2; char* s;
  if (! ifarg(1)) {printf("setdvi(v1,v2[,v3,flag]): v1:cell#s; v2:delays; v3:distal synapses\n"); return 0; }
  ip=IDP; pg=ip->pg; 
  if (ip->dead) return 0;
  dvt=vector_arg_px(1, &y);
  i=vector_arg_px(2, &d);
  s=ifarg(3)?(char*)calloc((j=vector_arg_px(3,&ds)),sizeof(char)):0x0;
  if(s) for(k=0;k<j;k++) s[k]=(char)ds[k];
  if (ifarg(4)) flag=(int)*getarg(4); else flag=0;
  if (i!=dvt || i==0 || (j>0 && j!=i)) {printf("setdvi() ERR vec sizes: %d %d %d\n",dvt,i,j); hxe();}
  w1=w2=0x0;
  if(ifarg(5) && wsetting!=1){printf("setdvi ERR: only use weight vecs when wsetting==1!\n"); hxe();}
  if(ifarg(5) && dvt!=vector_arg_px(5,&w1)){printf("setdvi ERR: wrong size for w1 vector!\n"); hxe();}
  if(ifarg(6) && dvt!=vector_arg_px(6,&w2)){printf("setdvi ERR: wrong size for w2 vector!\n"); hxe();}
  setdvi2(y,d,s,dvt,flag,w1,w2);
  }
ENDVERBATIM
}

VERBATIM


static int setdvi2 (double *y,double *d,char* s,int dvt,int flag,double* w1,double* w2) {
  int i,j,ddvi; double *db, *dbs, *w1s, *w2s; unsigned char pdead; unsigned int b,e; char* syns;
  Object *lb; Point_process *pnnt, **da, **das;
  ddvi=(int)DEAD_DIV;
  ip=IDP; pg=ip->pg;
  if(wsetting==1 && (!w1 || !w2)) {
    printf("setdvi2 ERR: wsetting==1, must provide w1,w2 arrays!\n");
    hxe();
  }
  if (flag==0) { b=0; e=dvt; 
    if (ip->dvi) { 
      free(ip->dvi); free(ip->del); free(ip->sprob); free(ip->syns); 
      ip->dvt=0; ip->dvi=(Point_process**)0x0; ip->del=(double*)0x0; ip->sprob=(char *)0x0; ip->syns=(char*)0x0;
      if(ip->wgain){free(ip->wgain); ip->wgain=0x0;}
      if(ip->peconv){free(ip->peconv); ip->peconv=0x0;}
      if(ip->piconv){free(ip->piconv); ip->piconv=0x0;}
      if(ip->pplasttau){free(ip->pplasttau);ip->pplasttau=0x0;}
      if(ip->pplastinc){free(ip->pplastinc);ip->pplastinc=0x0;}
      if(ip->pplastmaxw){free(ip->pplastmaxw);ip->pplastmaxw=0x0;}
      if(ip->pdope){free(ip->pdope);ip->pdope=0x0;}
      if(wsetting==1) freesywv(ip);
    } 
  } else { 
    if (ip->dvt==0) {
      ip->dvi=(Point_process**)0x0; ip->del=(double*)0x0; ip->sprob=(char *)0x0; ip->syns=(char*)0x0;
      ip->wgain=0x0; ip->peconv=0x0; ip->piconv=0x0;
      ip->pplasttau=0x0; ip->pplastinc=0x0; ip->pplastmaxw=0x0; ip->pdope=0x0;
      if(wsetting==1) freesywv(ip);
    }
    b=ip->dvt; 
    e=ip->dvt+dvt; 
  }
  da=(Point_process **)realloc((void*)ip->dvi,(size_t)(e*sizeof(Point_process *)));
  db=(double*)realloc((void*)ip->del,(size_t)(e*sizeof(double)));
  syns=(char*)realloc((void*)ip->syns,(size_t)(e*sizeof(char)));  
  if(wsetting==1) {
    w1s=(double*)realloc((void*)ip->syw1,(size_t)(e*sizeof(double)));
    w2s=(double*)realloc((void*)ip->syw2,(size_t)(e*sizeof(double)));
  }
  for (i=b,j=0;j<dvt;j++) { 
    
    if (!(lb=ivoc_list_item(pg->ce,(unsigned int)y[j]))) {
      printf("INTF6:callback %g exceeds %d for list ce\n",y[j],pg->cesz); hxe(); }
      pnnt=(Point_process *)lb->u.this_pointer;
      if (ddvi==1 || !(pdead=(*(id0**)&(pnnt->_prop->dparam[2]))->dead)) {
        da[i]=pnnt; db[i]=d[j]; syns[i]=s?s[j]:0; 
        if(wsetting==1){w1s[i]=w1[j]; w2s[i]=w2[j];}
        i++;
      }
  }
  if ((dvt=i)<e) { 
    da=(Point_process **)realloc((void*)da,(size_t)(dvt*sizeof(Point_process *)));
    db=(double*)realloc((void*)db,(size_t)(dvt*sizeof(double)));
    syns=(char*)realloc((void*)syns,(size_t)(dvt*sizeof(char)));
    if(wsetting==1) {
      w1s=(double*)realloc((void*)w1s,(size_t)(dvt*sizeof(double)));
      w2s=(double*)realloc((void*)w2s,(size_t)(dvt*sizeof(double)));
    }
  }
  ip->dvt=dvt; ip->del=db; ip->dvi=da; ip->syns=syns;
  if(wsetting==1){ip->syw1=w1s; ip->syw2=w2s;}
  if (flag!=1) finishdvi2(ip); 
}
ENDVERBATIM

VERBATIM


static int setdvi3 (double *y, double *d, char* s, int dvt, double* w1, double* w2) {
  int i,j,ddvi; double *db, *dbs, *w1s, *w2s; unsigned char pdead; unsigned int b,e; char* syns;
  Object *lb; Point_process *pnnt, **da, **das;
  ddvi=(int)DEAD_DIV;
  ip=qp; pg=ip->pg;
  e=dvt; 
  da=(Point_process **)realloc((void*)ip->dvi,(size_t)(e*sizeof(Point_process *)));
  db=(double*)realloc((void*)ip->del,(size_t)(e*sizeof(double)));
  syns=(char*)realloc((void*)ip->syns,(size_t)(e*sizeof(char)));  
  w1s=(double*)realloc((void*)ip->syw1,(size_t)(e*sizeof(double)));
  w2s=(double*)realloc((void*)ip->syw2,(size_t)(e*sizeof(double)));
  for (i=0,j=0;j<dvt;i++,j++) { 
    
    if (!(lb=ivoc_list_item(pg->ce,(unsigned int)y[j]))) {
      printf("INTF6:callback %g exceeds %d for list ce\n",y[j],pg->cesz); hxe(); }
      pnnt=(Point_process *)lb->u.this_pointer;
      if (ddvi==1 || !(pdead=(*(id0**)&(pnnt->_prop->dparam[2]))->dead)) {
        da[i]=pnnt; db[i]=d[j]; syns[i]=s?s[j]:0; 
        w1s[i]=w1[j]; w2s[i]=w2[j];
      }
  }
  ip->dvt=dvt; ip->del=db; ip->dvi=da; ip->syns=syns;
  ip->syw1=w1s; ip->syw2=w2s;
  finishdvi2(ip); 
}
ENDVERBATIM



PROCEDURE prune () {
  VERBATIM 
  {
  id0* ppost; double *x, p; int nx,j,potype;
  ip=IDP; pg=ip->pg;
  if (hoc_is_double_arg(1)) { 
    p=*getarg(1);
    if (p<0 || p>1) {printf("INTF6:pruneERR0:need # [0,1] to prune [ALL,NONE]: %g\n",p); hxe();}
    if (p==1.) printf("INTF6pruneWARNING: pruning 100% of cell %d\n",ip->id);
    if (verbose && ip->dvt>dscrsz) {
      printf("INTF6pruneB:Div exceeds dscrsz: %d>%d\n",ip->dvt,dscrsz); hxe(); }
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
        ppost=*((id0**) &((ip->dvi[j]->_prop->dparam)[2])); 
        if (ppost->type==potype && dscr[j]<p) ip->sprob[j]=0; 
      }
    }
  } else { 
    if (verbose) printf("INTF6 WARNING prune(vec) deprecated: use intf.sprob(vec) instead\n");
    nx=vector_arg_px(1,&x);
    if (nx!=ip->dvt) {printf("INTF6:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
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
  if (nx!=ip->dvt) {printf("INTF6:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
  if (ifarg(2)) { 
    if (!hoc_is_str_arg(2)) { printf("INTF6 sprob()ERRA: only legit 2nd arg is 'GET'\n"); hxe();
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
  ip=IDP; pg=ip->pg;
  nx=vector_arg_px(1,&x);
  ny=vector_arg_px(2,&y);
  if (ifarg(3)) off=(unsigned char)*getarg(3); else off=0;
  for (i=0;i<nx;i++) { 
    lop(pg->ce,(unsigned int)x[i]); 
    dvt=qp->dvt; das=qp->dvi;
    for (j=0;j<dvt;j++) {
      ip=*((id0**) &((das[j]->_prop->dparam)[2])); 
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

int gsort2 (double *db, Point_process **da,int dvt,double *dbs, Point_process **das) {
  int i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
  }
}

int gsort3 (double *db, Point_process **da,char* syns,int dvt,double *dbs, Point_process **das,char* synss) {
  int i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
    synss[i]=syns[scr[i]];
  }
}


int gsort5 (double *db, Point_process **da, char* syns, double* w1,double* w2, int dvt,
           double *dbs, Point_process **das,char* synss,double* w1s,double* w2s) {
  int i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
    synss[i]=syns[scr[i]];
    w1s[i]=w1[scr[i]];
    w2s[i]=w2[scr[i]];
  }
}

static int freedvi2 (struct ID0* jp) {
  if (jp->dvi) {
    free(jp->dvi); free(jp->del); free(jp->sprob); free(jp->syns);
    if(jp->wgain){free(jp->wgain); jp->wgain=0x0;}
    if(jp->peconv){free(jp->peconv); jp->peconv=0x0;}
    if(jp->piconv){free(jp->piconv); jp->piconv=0x0;}
    if(ip->pplasttau){free(ip->pplasttau);ip->pplasttau=0x0;}
    if(ip->pplastinc){free(ip->pplastinc);ip->pplastinc=0x0;}
    if(ip->pplastmaxw){free(ip->pplastmaxw);ip->pplastmaxw=0x0;}
    if(ip->pdope){free(ip->pdope);ip->pdope=0x0;}
    jp->dvt=0; jp->dvi=(Point_process**)0x0; jp->del=(double*)0x0; jp->sprob=(char *)0x0; jp->syns=(char *)0x0;
  }
}
ENDVERBATIM

PROCEDURE freedvi () {
  VERBATIM
  { 
    id0 *jp;
    jp=IDP;
    freedvi2(jp);
  }
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


FUNCTION mywmat () {
  VERBATIM {
  int i,j,k;
  i=(int)*getarg(1);
  if(i<0 || i>=CTYPi){printf("mywmat ERR: arg 1=%d out of bounds (0,%d]\n",i,CTYPi); return -1;}
  j = (int)*getarg(2);
  if(j<0 || j>=CTYPi){printf("mywmat ERR: arg 2=%d out of bounds (0,%d]\n",j,CTYPi); return -1;}
  k = (int)*getarg(3);
  if(k<0 || k>=STYPi){printf("mywmat ERR: arg3=%d out of bounds (0,%d]\n",k,STYPi); return -1;}
  return WMAT(i,j,k);
  }
  ENDVERBATIM  
}


PROCEDURE mywmatpr () {
  VERBATIM {
  double wm;
  int i,j,k;
  char *ct1,*ct2;
  ip=IDP; pg=ip->pg;
  for(i=0;i<CTYPi;i++) if(ctt(i,&ct1)!=0) {
    for(j=0;j<CTYPi;j++) if(ctt(j,&ct2)!=0) {
      for(k=0;k<STYPi;k++) {
        if((wm=WMAT(i,j,k))>0) {
          printf("wmat[%s][%s][%d]=%g\n",ct1,ct2,k,wm);
        }
      }      
    }
  }
  }
  ENDVERBATIM
}



PROCEDURE cinit () {
  VERBATIM {
  Symbol *sym; int i,j; unsigned int sz,colid; char *name;

  pg=(postgrp *)calloc(1,sizeof(postgrp));
  colid = (int)*getarg(2);

  if(ppg==0x0) { 
    ippgbufsz = 5;
    ppg = (postgrp**) calloc(ippgbufsz,sizeof(postgrp*));
    inumcols = 1;
  } else inumcols++;

  if(inumcols >= ippgbufsz) { 
    ippgbufsz *= 2;
    ppg = realloc((void*)ppg,(size_t)ippgbufsz*sizeof(postgrp*));
  }
  ppg[inumcols-1] = pg;
  pg->col = colid;
  pg->ce =  *hoc_objgetarg(1);

  sym = hoc_lookup("CTYP"); 
  CTYP = (*(hoc_objectdata[sym->u.oboff].pobj));

  if (installed==2.0) { 
    sz=ivoc_list_count(pg->ce);
    if (sz==pg->cesz && colid==0) printf("\t**** INTF6 WARNING cesz unchanged: INTF6(s) created off-list ****\n");
  } else installed=2.0;
  pg->cesz = ivoc_list_count(pg->ce); if(verbose) printf("cesz=%d\n",pg->cesz);
  pg->lastspk = calloc(pg->cesz,sizeof(double)); 
  
  CTYPi=HVAL("CTYPi"); STYPi=HVAL("STYPi"); dscrsz=HVAL("scrsz"); dscr=HPTR("scr");
  
  pg->ix = hoc_pgetarg(3);
  pg->ixe = hoc_pgetarg(4);
  pg->numc = hoc_pgetarg(5); 
  if(verbose){printf("CTYPi=%d\n",CTYPi);
    for(i=0;i<CTYPi;i++) printf("ix[%d]=%g, ixe[%d]=%g\n",i,pg->ix[i],i,pg->ixe[i]);}
  if (!pg->ce) {printf("INTF6 cinit() ERRA: ce not found\n"); hxe();}
  if (ivoc_list_count(CTYP)!=CTYPi){
    printf("INTF6 cinit() ERRB: %d %d\n",ivoc_list_count(CTYP),CTYPi); hxe(); }
  for (i=0;i<pg->cesz;i++) { lop(pg->ce,i); qp->pg=pg; } 
  printf("Checking for possible seg error in double arrays: CTYPi==%d: ",CTYPi);
  printf("%d %g\n",dscrsz,dscr[dscrsz-1]); 
  for (i=0,j=0;i<CTYPi;i++) if (ctt(i,&name)!=0) {
    cty[j]=i; CNAME[j]=name; ctymap[i]=j;
    j++;
    if (j>=CTYPp) {printf("jitcondiv() INTERRA\n"); hxe();}
  }
  CTYN=j; 
  for (i=0;i<CTYN;i++) printf("%s(%d)=%g ",CNAME[i],cty[i],NUMC(cty[i]));
  printf("\n%d cell types being used in col %d\n",CTYN,colid);
  }
  ENDVERBATIM  
}


PROCEDURE jitcondiv () {
  VERBATIM {
  Symbol *sym; int i,j; unsigned int sz,colid; char *name;

  pg=(postgrp *)calloc(1,sizeof(postgrp));
  colid = (int)*getarg(2);

  if(ppg==0x0) { 
    ippgbufsz = 5;
    ppg = (postgrp**) calloc(ippgbufsz,sizeof(postgrp*));
    inumcols = 1;
  } else inumcols++;

  if(inumcols >= ippgbufsz) { 
    ippgbufsz *= 2;
    ppg = realloc((void*)ppg,(size_t)ippgbufsz*sizeof(postgrp*));
  }
  ppg[inumcols-1] = pg;
  pg->col = colid;
  pg->ce =  *hoc_objgetarg(1);

  sym = hoc_lookup("CTYP"); 
  CTYP = (*(hoc_objectdata[sym->u.oboff].pobj));

  if (installed==2.0) { 
    sz=ivoc_list_count(pg->ce);
    if (sz==pg->cesz && colid==0) printf("\t**** INTF6 WARNING cesz unchanged: INTF6(s) created off-list ****\n");
  } else installed=2.0;
  pg->cesz = ivoc_list_count(pg->ce); if(verbose) printf("cesz=%d\n",pg->cesz);
  pg->lastspk = calloc(pg->cesz,sizeof(double)); 

  
  CTYPi=HVAL("CTYPi"); STYPi=HVAL("STYPi"); dscrsz=HVAL("scrsz"); dscr=HPTR("scr");

  
  pg->ix = hoc_pgetarg(3);
  pg->ixe = hoc_pgetarg(4);

  if(verbose){printf("CTYPi=%d\n",CTYPi);
    for(i=0;i<CTYPi;i++) printf("ix[%d]=%g, ixe[%d]=%g\n",i,pg->ix[i],i,pg->ixe[i]);}

  pg->dvg = hoc_pgetarg(5); 
  pg->numc = hoc_pgetarg(6); 
  pg->wmat = hoc_pgetarg(7); 
  pg->wd0 = hoc_pgetarg(8); 
  pg->delm = hoc_pgetarg(9); 
  pg->deld = hoc_pgetarg(10); 

  if (!pg->ce) {printf("INTF6 jitcondiv ERRA: ce not found\n"); hxe();}
  if (ivoc_list_count(CTYP)!=CTYPi){
    printf("INTF6 jitcondiv ERRB: %d %d\n",ivoc_list_count(CTYP),CTYPi); hxe(); }
  for (i=0;i<pg->cesz;i++) { lop(pg->ce,i); qp->pg=pg; } 
  
  printf("Checking for possible seg error in double arrays: CTYPi==%d: ",CTYPi);
  
  printf("%d %d %d ",DVG(CTYPi-1,CTYPi-1),(int)pg->ix[CTYPi-1],(int)pg->ixe[CTYPi-1]);
  printf("%g %g ",WMAT(CTYPi-1,CTYPi-1,STYPi-1),WD0(CTYPi-1,CTYPi-1,STYPi-1));
  printf("%g %g ",DELM(CTYPi-1,CTYPi-1),DELD(CTYPi-1,CTYPi-1));
  printf("%d %g\n",dscrsz,dscr[dscrsz-1]); 
  for (i=0,j=0;i<CTYPi;i++) if (ctt(i,&name)!=0) {
    cty[j]=i; CNAME[j]=name; ctymap[i]=j;
    j++;
    if (j>=CTYPp) {printf("jitcondiv() INTERRA\n"); hxe();}
  }
  CTYN=j; 
  for (i=0;i<CTYN;i++) printf("%s(%d)=%g ",CNAME[i],cty[i],NUMC(cty[i]));
  printf("\n%d cell types being used in col %d\n",CTYN,colid);
  }
  ENDVERBATIM  
}


PROCEDURE jitrec () {
  VERBATIM {
  int i;
  ip=IDP; pg=ip->pg;
  if(verbose>1) printf("jitrec from col %d, ip=%p, pg=%p\n",ip->col,ip,pg);
  if (! ifarg(2)) { 
    pg->jrmax=0; pg->jridv=0x0; pg->jrtvv=0x0;
    return 0;
  }
  i =   vector_arg_px(1, &pg->jrid); 
  pg->jrmax=vector_arg_px(2, &pg->jrtv);
  pg->jridv=vector_arg(1); pg->jrtvv=vector_arg(2);
  pg->jrmax=vector_buffer_size(pg->jridv);
  if (pg->jrmax!=vector_buffer_size(pg->jrtvv)) {
    printf("jitrec() ERRA: not same size: %d %d\n",i,pg->jrmax); pg->jrmax=0; hxe(); }
  pg->jri=pg->jrj=0; 
  }
  ENDVERBATIM
}






FUNCTION scsv () {
  VERBATIM {
  int ty=4; int i,j; unsigned int cnt=0;
  ip=IDP; pg=ip->pg;
  name = gargstr(1);
  if ( !(wf1 = fopen(name,"w"))) { printf("Can't open %s\n",name); hxe(); }
  fwrite(&pg->cesz,sizeof(int),1,wf1);
  fwrite(&ty,sizeof(int),1,wf1);
  for (i=0,j=0;i<pg->cesz;i++,j++) { 
    lop(pg->ce,i); 
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
  ip=IDP; pg=ip->pg;
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
  if (ny) max=ny; else if (max==0) max=pg->cesz; else max+=1; 
  if (nx && nx!=max-min) {
    printf("INTF6 spkcnt() ERR: Vectors not same size %d %d\n",nx,max-min);hxe();}
  for  (i=min, sum=0;i<max;i++) { 
    if (ny) lop(pg->ce,(int)y[i]); else lop(pg->ce,i);
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
    ip=IDP; pg=ip->pg;
    for (i=1;i<=3;i++) a[i]=(int)*getarg(i);
    printf("CTYPi: %d, STYPi: %d, ",CTYPi,STYPi);
    
    printf("wmat: %g, wd0: %g\n",WMAT(a[1],a[2],a[3]),WD0(a[1],a[2],a[3]));
  }
  ENDVERBATIM  
}


PROCEDURE randspk () {
  VERBATIM 
  ip=IDP; pg=ip->pg;
  if (ip->rvi > ip->rve) { 
    ip->input=0;           
    
                
                
  } else if (t==0) {     
    nxt=pg->vsp[ip->rvi];
    EXSY=pg->sysp[ip->rvi]; 
    WEX=pg->wsp[ip->rvi++]; 
  } else {     
    while ((nxt=pg->vsp[ip->rvi++]-t)<=1e-6) { 
      if (ip->rvi-1 > ip->rve) { printf("randspk() ERRA: "); chk(2.); hxe(); }
    }
    EXSY=pg->sysp[ip->rvi-1]; 
    WEX=pg->wsp[ip->rvi-1]; 
  }
  ENDVERBATIM
  
}


PROCEDURE vers () {
  printf("$Id
}


VERBATIM
double val (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);

  vii[5]=AHP*EXP(-(xx - ta)/tauahp);
  vii[8]=VAM2*EXP(-(xx -ta)/tauAM2);
  vii[9]=VNM2*EXP(-(xx - ta)/tauNM2);
  vii[10]=VGA2*EXP(-(xx - ta)/tauGA2);
  vii[6]=vii[1]+vii[2]+vii[3]+vii[4]+vii[5]+vii[8]+vii[9]+vii[10];
  vii[7]=VTH + (VTHR-VTH)*EXP(-(xx-trrs)/tauRR) - RMP; 
}
ENDVERBATIM


VERBATIM
double valps (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);

  vii[8]=VAM2*EXP(-(xx - ta)/tauAM2);
  vii[9]=VNM2*EXP(-(xx - ta)/tauNM2);
  vii[10]=VGA2*EXP(-(xx - ta)/tauGA2);
  vii[6]=vii[1]+vii[2]-vii[3]+vii[8]+vii[9]-vii[10];
}
ENDVERBATIM


PROCEDURE record () {
  VERBATIM {
  int i,j,k,nz; double ti;
  vp = SOP;
  if(!vp) {printf("**** record ERRA: vp=NULL!\n"); return 0;}
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
  VERBATIM { 
  vp = SOP;
  record();
  if (vp->p > vp->size || vp->vvo[6]==0) return 0; 
  if (vp->p > 0) {
    if (vp->vvo[0]!=0x0) vp->vvo[0][vp->p-1]=_lx;
    vp->vvo[6][vp->p-1]=spkht; 
  } else {
    if (vp->vvo[0]!=0x0) vp->vvo[0][0]=_lx; 
    vp->vvo[6][0]=spkht; 
  }
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
    } else printf("INTF6 recclr ERR: nil pointer\n");
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
  } else printf("INTF6 recfree ERR: nil pointer\n");
  IDP->record=0;
  ENDVERBATIM
}






PROCEDURE initvspks () {
  VERBATIM
  {int max, i,err;
    double last,lstt;
    ip=IDP; pg=ip->pg;
    if (! ifarg(1)) {printf("Return initvspks(indices,times,weights,syntypes)\n"); return 0.;}
    if(verbose>1) printf("initvspks: col=%d, ip=%p, pg=%p, pg->isp=%p\n",ip->col,ip,pg,pg->isp);
    if (pg->isp!=NULL) clrvspks();
    ip=IDP; pg=ip->pg; err=0;
    i = vector_arg_px(1, &pg->isp); 
    max=vector_arg_px(2, &pg->vsp);
    if (max!=i) {err=1; printf("initvspks ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("initvspks ERR: vec not initialized\n");}
    max=vector_arg_px(3, &pg->wsp);
    if (max!=i) {err=1; printf("initvspks ERR: 3rd vec is of different size\n");}
    max=vector_arg_px(4, &pg->sysp);
    if (max!=i) {err=1; printf("initvspks ERR: 4th vec is of different size\n");}
    pg->vspn=max;
    if (!pg->ce) {printf("Need global ce for initvspks() since intf.mod501\n"); hxe();}
    for (i=0,last=-1; i<max; ) { 
      if (pg->isp[i]!=last) { 
        lop(pg->ce,(unsigned int)pg->isp[i]);
        qp->rvb=qp->rvi=i;
        qp->vinflg=qp->input=1;
        last=pg->isp[i];
        lstt=pg->vsp[i];
        i++;
      }
      for (; i<max && pg->isp[i] == last; i++) { 
        if (pg->vsp[i]<=lstt) { pg->vsp[i]=lstt+0.00001; 
          printf("initvspks ERR: nonmonotonic for cell#%d: %g %g\n",qp->id,lstt,pg->vsp[i]); }
          lstt=pg->vsp[i];
      }
      qp->rve=i-1;
      if (subsvint>0) { 
        pg->vsp[qp->rve] = pg->vsp[qp->rvb]+subsvint;
        pg->wsp[qp->rve] = pg->wsp[qp->rvb];
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
    printf("WARNING: This routine appears to be defunct -- please check code in intf6.mod\n");
    if (! ifarg(1)) {printf("Return shock(ivspks,vspks,wvspks)\n"); return 0.;}
    ip=IDP; pg=ip->pg; err=0;
    i = vector_arg_px(1, &isp); 
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("shock ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("shock ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("shock ERR: 3rd vec is of different size\n");}
    pg->vspn=max;
    if (!pg->ce) {printf("Need global ce for shock()\n"); hxe();}
    for (i=0,last=-1; i<max; ) { 
      if (isp[i]!=last) { 
        lop(pg->ce,(unsigned int)isp[i]);
        WEX=-1e9; 
        EXSY=AM;  
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
 ip=IDP; pg=ip->pg;
 if(verbose>1) printf("clrvspks: col=%d, ip=%p, pg=%p, pg->isp=%p\n",ip->col,ip,pg,pg->isp);
 for (i=0; i<pg->cesz; i++) {
   lop(pg->ce,i);
   qp->vinflg=0;
 }   
 }
 ENDVERBATIM
}




PROCEDURE trvsp ()
{
  VERBATIM 
  int i, flag; 
  double ind, t0;
  ip=IDP; pg=ip->pg;
  flag=(int) *getarg(1);
  if (subsvint==0.) {printf("trvsp"); return(0.);}
  ind=pg->isp[0]; t0=pg->vsp[0];
  if (flag==1) {
    for (i=0; i<pg->vspn; i++) {
      if (pg->isp[i]!=ind) {
        pg->vsp[i-1]=1.e9;
        ind=pg->isp[i];
      }
    }
    pg->vsp[pg->vspn-1]=1.e9;
  } else if (flag==2) {
    for (i=0; i<pg->vspn; i++) {
      if (pg->isp[i]!=ind) {
        pg->vsp[i-1]=t0+subsvint;
        ind=pg->isp[i]; t0=pg->vsp[i];
      }
    }
    pg->vsp[pg->vspn-1]=t0+subsvint;
  } else {printf("trvsp flag %d not recognized\n",flag); hxe();}
  ENDVERBATIM
}




PROCEDURE initjttr () {
  VERBATIM 
  {int max, i, err=0;
    ip=IDP; pg=ip->pg;
    pg->jtpt=0;
    if (! ifarg(1)) {printf("Return initjttr(vec)\n"); return(0.);}
    max=vector_arg_px(1, &jsp);
    if (max==0) {err=1; printf("initjttr ERR: vec not initialized\n");}
    for (i=0; i<max; i++) if (jsp[i]<=0) {err=1;
      printf("initjttr ERR: vec should be >0: %g\n",jsp[i]);}
    if (err) { jsp=nil; pg->jtmax=0.; return(0.); }
    if (max != pg->jtmax) {
      printf("WARNING: resetting jtmax_INTF6 to %d\n",max); pg->jtmax=max; }
  }
  ENDVERBATIM
}


VERBATIM



static id0* getlp (Object *ob, unsigned int i) {
  Object *lb; id0* myp;
  lb = ivoc_list_item(ob, i);
  if (! lb) { printf("INTF6:getlp %d exceeds %d for list ce\n",i,pg->cesz); hxe();}
  pmt=ob2pntproc(lb);
  myp=*((id0**) &((pmt->_prop->dparam)[2])); 
  return myp;
}



static id0* lop (Object *ob, unsigned int i) {
  Object *lb;
  lb = ivoc_list_item(ob, i);
  if (! lb) { printf("INTF6:lop %d exceeds %d for list ce\n",i,pg->cesz); hxe();}
  pmt=ob2pntproc(lb);
  qp=*((id0**) &((pmt->_prop->dparam)[2])); 
  return qp;
}


static id0* lopr (Object *ob, unsigned int i) {
  id0* myp;
  myp = lop(ob,i);
  _hoc_setdata(pmt); 
  return myp;
}


int stoppo () {
}


static int ctt (unsigned int i, char** name) {
  Object *lb;
  if (NUMC(i)==0) return 0; 
  lb = ivoc_list_item(CTYP, i);
  if (! lb) { printf("INTF6:ctt %d exceeds %d for list CTYP\n",i,CTYPi); hxe();}
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
  Object *ob; int num,i,ii,j,k,si,nx;  double *vvo[7], *par; void *vv[7];
  ob = *(hoc_objgetarg(1));
  si=(int)*getarg(2);
  num = ivoc_list_count(ob);
  if (num!=7) { printf("INTF6 lof ERR %d>7\n",num); hxe(); }
  for (i=0;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);
    if (i==0) nx=j;
    if (j!=nx) { printf("INTF6 lof ERR %d %d\n",j,nx); hxe(); }
  }
  
  
  
  
  
  
  
 }
ENDVERBATIM
}



PROCEDURE initinvl () {
  printf("initinvl() NOT BEING USED\n")
}


FUNCTION invlflag () {
  VERBATIM
  ip=IDP; pg=ip->pg;
  if (ip->invl0==1 && invlp==nil) { 
    printf("INTF6 invlflag ERR: pointer not initialized\n"); hoc_execerror("",0); 
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
    printf("INTF6 record ERR: pointer not initialized\n"); hoc_execerror("",0); 
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
  ip=IDP; pg=ip->pg;
  printf("ID:%d; typ: %d; rec:%d wrec:%d inp:%d jtt:%d invl:%d\n",ip->id,ip->type,ip->record,ip->wrec,ip->input,ip->jttr,ip->invl0);
  if (lfg==1) {
    if (SOP!=nil) {
      vp = SOP;
      printf("p %d size %d tg %g\n",vp->p,vp->size,tg);
      for (i=0;i<NSV;i++) if (vp->vv[i]) printf("%d %x %x;",i,(unsigned int)vp->vv[i],(unsigned int)vp->vvo[i]);
    } else printf("Recording pointers not initialized");
  }
  if (lfg==2) { 
    printf("Global vectors for input and jitter (jttr): \n");
    if (pg->vsp!=nil) printf("VSP: %x (%d/%d-%d)\n",(unsigned int)pg->vsp,ip->rvi,ip->rvb,ip->rve); else printf("no VSP\n");
    if (jsp!=nil) printf("JSP: %x (%d/%d)\n",(unsigned int)jsp,pg->jtpt,pg->jtmax); else printf("no JSP\n");
  }
  if (lfg==3) { 
    if (pg->vsp!=nil) { printf("VSP: (%d/%d-%d)\n",ip->rvi,ip->rvb,ip->rve); 
      for (i=ip->rvb;i<=ip->rve;i++) printf("%d:%g  ",i,pg->vsp[i]);
      printf("\n");
    } else printf("no VSP\n");
  }
  if (lfg==4) {  
  }
  if (lfg==5) { 
    printf("wwpt %d wwsz %d\n WW vecs: ",wwpt,wwsz);
    printf("wwwid %g wwht %d nsw %g\n WW vecs: ",wwwid,(int)wwht,nsw);
    for (i=0;i<NSW;i++) printf("%d %x %x;",i,(unsigned int)ww[i],(unsigned int)wwo[i]);
  }}
  ENDVERBATIM
}


FUNCTION pid () {
  VERBATIM 
  printf("INTF6%d(%d/%d@%g) ",IDP->id,IDP->type,IDP->col,t);
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


FUNCTION gid () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->gid = (unsigned int) *getarg(1);
  _lgid = (double)ip->gid;
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
  if (i==-1) {printf("INTF6 record ERR %s not recognized\n",name); hoc_execerror("",0); }
  vp->vv[i]=vector_arg(2);
  vector_arg_px(2, &(vp->vvo[i]));
  if (vp->size==0) { vp->size=(unsigned int)vector_buffer_size(vp->vv[i]);
  } else if (vp->size != (unsigned int)vector_buffer_size(vp->vv[i])) {
    printf("INTF6 initrec ERR vectors not all same size: %d vs %d",vp->size,vector_buffer_size(vp->vv[i]));
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
  } else if (strcmp(name,"AHP")==0) { _li=5.;
  } else if (strcmp(name,"V")==0) { _li=6.;
  } else if (strcmp(name,"VM")==0) { _li=6.; 
  } else if (strcmp(name,"VTHC")==0) { _li=7.;
  } else if (strcmp(name,"VAM2")==0) { _li=8.;
  } else if (strcmp(name,"VNM2")==0) { _li=9.;
  } else if (strcmp(name,"VGA2")==0) { _li=10.;
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
  else if (i==5) printf("AHP\n");
  else if (i==6) printf("V\n");
  else if (i==7) printf("VTHC\n");
  else if (i==8) printf("VAM2\n");
  else if (i==9) printf("VNM2\n");
  else if (i==10) printf("VGA2\n");
  ENDVERBATIM
}


PROCEDURE initwrec () {
  VERBATIM 
  {int i, k, num, cap;  Object* ob;
    ob =   *hoc_objgetarg(1); 
    num = ivoc_list_count(ob);
    if (num>NSW) { printf("INTF6 initwrec() WARN: can only store %d ww vecs\n",NSW); hxe();}
    nsw=(double)num;
    for (k=0;k<num;k++) {
      cap = list_vector_px2(ob, k, &wwo[k], &ww[k]);
      if (k==0) wwsz=cap; else if (wwsz!=cap) {
        printf("INTF6 initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,cap); hxe(); }
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
        max=(int)_tmin_popspk; 
        k=-(int)floor((_lte-rebeg)/vdt+0.5);
        for (j= -max;j<=max && k+j>0 && k+j<wwsz;j++) {
          wwo[wrp][k+j] += scale*_t_Psk[j+max]; 
        }
      }
    } else if (twg>=t) { return 0;
    } else {
      for (ti=twg,k=(int)floor((twg-rebeg)/vdt+0.5);ti<=t && k<wwsz;ti+=vdt,k++) { 
        valps(ti,twg);  
        wwo[wrp][k]+=vii[6]*lfpscale;
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
    } else {printf("INTF6 wrec ERR flag(0/1) %d\n",ip->wrec); hxe();
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
  ip=IDP; pg=ip->pg;
  if (pg->jtmax>0 && pg->jtpt>=pg->jtmax) {  
    pg->jtpt=0;
    printf("Warning, cycling through jttr vector at t=%g\n",t);
  }
  if (pg->jtmax>0) _ljttr = jsp[pg->jtpt++]; else _ljttr=0;
  ENDVERBATIM
}


PROCEDURE global_init () {
  popspk(0) 
  VERBATIM 
  { int i,j,k,c; double stt[3];
  if (nsw>0. && wwo[0]!=0) { 
    printf("Initializing ww to record for %g (%g)\n",vdt*wwsz,vdt);
    wwpt=0;
    for (k=0;k<(int)nsw;k++) {
      vector_resize(ww[k], wwsz);
      for (j=0;j<wwsz;j++) wwo[k][j]=0.;
    }
  }
  errflag=0;
  for (i=0;i<CTYN;i++) blockcnt[cty[i]]=spikes[cty[i]]=0;
  for(c=0;c<inumcols;c++) {
    pg=ppg[c]; if(!pg) continue;
    if (pg->jridv) { pg->jri=pg->jrj=0; vector_resize(pg->jridv, pg->jrmax); vector_resize(pg->jrtvv, pg->jrmax); }
    pg->spktot=0;
    pg->jtpt=0;
    pg->eventtot=0;
  }
  }
  ENDVERBATIM
}

PROCEDURE global_fini () {
  VERBATIM
  int c,k;
  for (k=0;k<(int)nsw;k++) vector_resize(ww[k], (int)floor(t/vdt+0.5));
  for(c=0;c<inumcols;c++) {
    pg=ppg[c]; if(!pg) continue;
    if (pg->jridv && pg->jrj<pg->jrmax) {
      vector_resize(pg->jridv, pg->jrj); 
      vector_resize(pg->jrtvv, pg->jrj);
    }
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
  ip=IDP; pg=ip->pg;
  if (ip->vinflg==0 && pg->vsp==nil) { 
  } else if (ip->vinflg==1 && ip->rve==-1) {
    printf("INTF6 vinflag ERR: pointer not initialized\n"); hoc_execerror("",0); 
  } else if (ip->rve >= 0) { 
    if (pg->vsp==nil) {
      printf("INTF6 vinflag ERR1: pointer not initialized\n"); hoc_execerror("",0); 
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
  ip=IDP; pg=ip->pg;
  if (FLAG==OK) { 
    FLAG=0.;
    if (stoprun) {slowset=0; return 0.0;}
    if (IDP->dbx==-1)printf("slowset fi:%d ix:%d ss:%g delt:%g t:%g\n",fi,ix,slowset,delt,t);
    if (t>slowset || ix>=pg->cesz) {  
      printf("Slow-setting of flag %d finished at %g: (%d,%g,%g)\n",fi,t,ix,delt,slowset); 
      slowset=0.; return 0.0;
    }
    if (ix<pg->cesz) {
      lop(pg->ce,ix);
      (&qp->type)[fi]=((fi>=iflneg)?(char)x[ix]:(unsigned char)x[ix]);
      ix++;
      #if defined(t)
      net_send((void**)0x0, wts,tpnt,t+delt,OK); 
      #else
      net_send((void**)0x0, wts,tpnt,delt,OK);
      #endif
    }
    return 0.0;
  }  
  if (slowset>0 && ifarg(3)) {
    printf("INTF6 flag() slowset ERR; attempted set during slowset: fi:%d ix:%d ss:%g delt:%g t:%g",\
           fi,ix,slowset,delt,t); 
    return 0.0;
  }
  ip = IDP; setfl=ifarg(3); 
  if (ifarg(4)) { slowset=*getarg(4); delt=slowset/pg->cesz; slowset+=t; } 
  sf = gargstr(1);
  for (fi=0;fi<iflnum && strncmp(sf, &iflags[fi*4], 3)!=0;fi++) ; 
  if (fi==iflnum) {printf("INTF6 ERR: %s not found as a flag (%s)\n",sf,iflags); hxe();}
  if (ifarg(2)) {
    if (hoc_is_double_arg(2)) {  
      val=(unsigned char)*getarg(2);
      if (slowset) { 
        printf("NOT IMPLEMENTED\n"); 
      } else if (setfl) { 
        for (ix=0;ix<pg->cesz;ix++) { lop(pg->ce,ix); (&qp->type)[fi]=val; }
      } else {  
        (&ip->type)[fi]=((fi>=iflneg)?(char)val:val);
      }
    } else {
      nx=vector_arg_px(2,&x);
      if (nx!=pg->cesz) {
        if (setfl) { printf("INTF6 flag ERR: vec sz mismatch: %d %d\n",nx,pg->cesz); hxe();
        } else       x=vector_newsize(vector_arg(2),pg->cesz);
      }
      if (setfl && slowset) { 
        ix=0;
        lop(pg->ce,ix);
        (&qp->type)[fi]=((fi>=iflneg)?(char)x[ix]:(unsigned char)x[ix]);
        ix++;
        #if defined(t)
        net_send((void**)0x0, wts,tpnt,t+delt,OK); 
        #else
        net_send((void**)0x0, wts,tpnt,delt,OK);
        #endif
      } else for (ix=0;ix<pg->cesz;ix++) { 
        lop(pg->ce,ix); 
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
  int i; double *x, sum; void *voi;
  ip = IDP; pg=ip->pg;
  voi=vector_arg(1);  x=vector_newsize(voi,pg->cesz);
  for (i=0,sum=0;i<pg->cesz;i++) { lopr(pg->ce,i); 
    x[i]=spck;
    sum+=spck;
  }
  _lallspck=sum;
  ENDVERBATIM
}


PROCEDURE resetall () {
  VERBATIM
  int ii,i; unsigned char val;
  ip=IDP; pg=ip->pg;
  if(verbose>1) printf("resetall: ip=%p, col=%d, pg=%p\n",ip,pg->col,pg);
  for (i=0;i<pg->cesz;i++) { lopr(pg->ce,i);
    Vm=RMP; VAM=0; VNM=0; VGA=0; AHP=0; invlt=-1; VAM2=0; VNM2=0; VGA2=0;
    t0=t; trrs=t; twg = t; cbur=0; spck=0; refractory=0; VTHC=VTHR=VTH; 
  }
  ENDVERBATIM
}



FUNCTION floc () {
  VERBATIM
  double x,y,z,r,min,rad, *ix, *dd, *tdy; int ii,i,n,cnt,ty,tvf; void *voi, *vod, *voty;
  cnt=0; n=1000; r=0; z=ty=1e9; tvf=0;
  ip = IDP; pg=ip->pg;
  x = *getarg(1);
  y = *getarg(2);
  i=3;
  if (ifarg(i)) if (hoc_is_double_arg(i)) { z=*getarg(3); i++; }
  if (ifarg(i)) {
    voi=vector_arg(i++); ix=vector_newsize(voi,n); 
    vod=vector_arg(i++); dd=vector_newsize(vod,n); 
    r= *getarg(i++); 
  }
  if (ifarg(i)) if (hoc_is_double_arg(i)) ty= *getarg(7); else { 
    tvf=1; voty=vector_arg(i++); tdy=vector_newsize(voty,n); 
  } 
  r*=r; 
  for (i=0,min=1e9,ii=-1;i<pg->cesz;i++) { qp=lopr(pg->ce,i); 
    if (ty!=1e9 && ((ty>=0 && ty!=qp->type) || (ty==-1 && qp->inhib==1) || (ty==-2 && qp->inhib==0))) continue;
    rad=(x-xloc)*(x-xloc)+(y-yloc)*(y-yloc)+(z==1e9?0.:((z-zloc)*(z-zloc))); 
    if (r>0 && rad<r) {
      
      if (cnt>=n) { 
        ix=vector_newsize(voi,n*=2); 
        dd=vector_newsize(vod,n);
        if (tvf) tdy=vector_newsize(voty,n);
      }
      ix[cnt]=(double)i;
      dd[cnt]=sqrt(rad); 
      if (tvf) tdy[cnt]=(double)qp->type;
      cnt++;
    } else if (rad<min) { min=rad; ii=i; }
  }
  if (r>0) { 
    ix=vector_newsize(voi,cnt); dd=vector_newsize(vod,cnt); 
    if (tvf) tdy=vector_newsize(voty,cnt);
    _lfloc=(double)cnt; } else {
    _lfloc=(double)ii;  } 
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

PROCEDURE ESINo (x) {
  TABLE ESIN FROM 0 TO 2*PI WITH 3000 
  ESIN = sin(x)
}

FUNCTION rates (vv) {
  
  rates = maxnmc / (1 + exp(0.062 (/mV) * -vv) * ( (mg / mg0) ) )
}









FUNCTION dopelearn () {
  VERBATIM
  Point_process *pnnt;
  int i , iCell, pot;
  double tmp,maxw,inc,d,tau,pdopet;
  if(seadsetting!=3.) return 0.0; 
  pot = (int) *getarg(1); 
  ip=IDP; pg=ip->pg;
  for(iCell=0;iCell<pg->cesz;iCell++){
    lop(pg->ce,iCell);
    if(!qp->dvt)continue; 
    for(i=0;i<qp->dvt;i++){
      pdopet = fabs(qp->pdope[i]); 
      d = t - pdopet; 
      if(qp->pdope[i] > -1e9 && d <= maxeligtrdur ) { 
        tmp = qp->wgain[i]; 
        maxw = qp->pplastmaxw[i]; 
        tau = qp->pplasttau[i]; 
        if( ! ( inc = qp->pplastinc[i] ) ) continue; 
        if(pot>0) { 
          if(qp->pdope[i] >= 0) { 
            if(SOFTSTDP) inc *= (1.0 - tmp / maxw); 
            if (EXPELIGTR) 
              qp->wgain[i] += EPOTW * inc * exp( -d / tau ); 
            else
              qp->wgain[i] += EPOTW * inc;
          } else { 
            if(SOFTSTDP) inc *= (tmp / maxw); 
            if (EXPELIGTR) 
              qp->wgain[i] -= EDEPW * inc * exp( -d / tau ); 
            else
              qp->wgain[i] -= EDEPW * inc;
          }
        } else { 
            if(qp->pdope[i] >= 0) { 
              if(SOFTSTDP) inc *= (tmp / maxw); 
              if (EXPELIGTR) 
                qp->wgain[i] -= EDEPW * inc * exp( -d / tau ); 
              else
                qp->wgain[i] -= EDEPW * inc;
            } else { 
              if(SOFTSTDP) inc *= (1.0 - tmp / maxw); 
              if (EXPELIGTR) 
                qp->wgain[i] += EPOTW * inc * exp( -d / tau ); 
              else
                qp->wgain[i] += EPOTW * inc;
            }
        }
        
        if(qp->wgain[i]<0.) qp->wgain[i]=0.; else if(!SOFTSTDP && qp->wgain[i]>maxw) qp->wgain[i]=maxw;
        if(reseteligtr) qp->pdope[i] = -1e9; 
      }
    }
  }
  return 1.0;
  ENDVERBATIM
}


PROCEDURE setdeletion () {
  VERBATIM
  
  
  
  
  
  
  
  double x = *getarg(1);
  if (x <= 0) {
    dynamicdel = 0; 
  } else {
    dynamicdel = 1; 
    delspeed = x; 
    printf("Set dynamic deletion rate constant = %e\n", delspeed);
  }
  ENDVERBATIM
}

VERBATIM
void dynamicdelete (double time) {
  
  
  
  
  
  mcell_ran4(&sead, dscr, 1, 1.0); 
  double p = dscr[0]; 

  double difference = ip->activity / ip->goal_activity; 

  
  
  
  
  
  
  
  double x = difference - 2.0;
  if (x < 0) {
    x = 0; 
  }
  double threshold =  x * x * delspeed * ((time - ip->lastupdate) / 1000.0);

  if (p < threshold) {
    printf("p = %e, threshold = %e from x^2 * delspeed * timegap:\nx = %e, x^2 = %e, delspeed = %e, x^2*delspeed = %e, timegap (s) = %e\n", p, threshold, x, x*x, delspeed, x*x*delspeed, (time-ip->lastupdate)/1000.0);

    ip->dead = 1; 
    printf("Cell %d has just died (scalefactor = %f)\n\n", ip->id, ip->scalefactor);
  }
}
ENDVERBATIM

VERBATIM
double get_avg_activity () {
  
  
  
  
  
  
  
  
  
  
  return ip->spkcnt / t;
}
ENDVERBATIM

VERBATIM 
void raise_activity_sensor (double time) {
  
  
  
  
  
  
  ip->activity = ip->activity + (-ip->activity + 1.0) / activitytau;

  
  
  
  
  
  
  
  
  
  
}
ENDVERBATIM

VERBATIM
void decay_activity_sensor (double time) {
  
  
  
  

  
  ip->activity = ip->activity * exp(-activityoneovertau * (time - ip->lastupdate));
}
ENDVERBATIM

VERBATIM
void update_scale_factor (double time) {
  

  
  double err = ip->goal_activity - ip->activity;

  
  
  
  
    
  
  
    
  

  
  ip->scalefactor += (activitybeta * ip->scalefactor * err + activitygamma * ip->scalefactor * ip->activity_integral_err);

  
  if (ip->scalefactor > ip->max_scale) {
    ip->scalefactor = ip->max_scale;
  }

  
  double timecorrection = time - ip->lastupdate;
  
  
  

  ip->activity_integral_err += (err * timecorrection);
  
  
  
  
}
ENDVERBATIM



FUNCTION scalefactor () {
  VERBATIM
  if (ifarg(1)) IDP->scalefactor = *getarg(1);
  return IDP->scalefactor; 
  ENDVERBATIM
}



FUNCTION activity () {
  VERBATIM
  if (ifarg(1)) IDP->activity = *getarg(1);
  return IDP->activity; 
  ENDVERBATIM
}



FUNCTION goalactivity () {
  VERBATIM
  if (ifarg(1)) IDP->goal_activity = *getarg(1);
  return IDP->goal_activity;   
  ENDVERBATIM
}


FUNCTION isdead() {
  VERBATIM
  return IDP->dead; 
  ENDVERBATIM
}