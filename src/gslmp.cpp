//2013.11.04 G.S. Jung
#include "gslmp.h"
using namespace LAMMPS_NS;
/*gsLMP::gsLMP(){
  lammps=NULL;
  thermo=NULL;

  lmpinput=NULL;
  cmd=NULL;
  fcc=NULL;
  graphene=NULL;
  cells=NULL;
  atominfo=NULL;

  metal="";
  potential="";
  vmdout="";
  cfgout="";
  timestep = 0.0;
  deformStep=0;
  thermoStep=0;
  rescaleStep=0;
  outputStep=0;
  scanStep=0;
  crackL=0.0;
  lattice=0.0;
  erate=0.0;
  mstep=0;
  lx0=0.0;
  ly0=0.0;
  lz0=0.0;

  on_crack=false;
  on_deform=false;
  on_strain=false;
  on_stmin=false;
  on_vmd=false;
  on_cfg=false;
  on_tslaw=false;

}

gsLMP::~gsLMP(){
  if(lammps!=NULL)delete lammps;
  lammps=NULL;
  thermo=NULL;
  if(fcc!=NULL)delete fcc;
  if(cells!=NULL)delete cells;
  if(cmd!=NULL)delete cmd;
  atominfo=NULL;
  lmpinput=NULL;
  }*/



void gsLMP::SetVMD(int _ostep, string _s){
  on_vmd=true;
  vmdStep=_ostep;
  vmdout=_s;
}

void gsLMP::SetCFG(int _ostep, string _s){
  on_cfg=true;
  cfgStep=_ostep;
  cfgout=_s;
}

void gsLMP::Equilibrium(int _tostep,double _temperature,double _pressure){
  ReadDat();
  if(on_tslaw){ cmd->Boundary("p s p");
    //cmd->Minimize();
  }
  else{ cmd->Boundary("p p p");
    //cmd->MinimizeISO();
  }

  cmd->TimeNeighbor(timestep,2.0);
  cmd->Thermos(thermoStep,"custom step temp pe press pxx pyy pzz lx ly lz");
  char com[200]={NULL};
  if(!on_tslaw) {

  }

  sprintf(com,"velocity all create %f 12345 mom yes rot no loop local",_temperature);
  lammps->input->one(com);
  
  char com2[200]={NULL};
  if(!on_tslaw) sprintf(com2,"all npt temp %f %f 1 iso %f %f 1 drag 1",_temperature,_temperature,_pressure,_pressure);
  else sprintf(com2,"all npt temp %f %f 1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure);
  
  cmd->Fix(com2);

  mstep = thermoStep;
  int tstep = _tostep/mstep+1;
  if(myrank==0) ofthermo.open(thermoout.c_str());

  cmd->Run(0);
  cmd->ResetTime();
  etotal0=0;
  for(int step=0;step<tstep;step++){//total iteration! 
    PrintThermo(step);
    cmd->Run(mstep);
  }

  if(myrank==0) ofthermo.close();

  cmd->UnFix(com2);  

}

void gsLMP::WriteRestart(string _readfile){
  char com[200]={NULL};
  sprintf(com,"write_restart %s",_readfile.c_str());
  lammps->input->one(com);
}

void gsLMP::Crack(int width,int length){
  fcc->RemoveCenterxc(width,length);  
  if(width>1)//fcc->RemoveCenterxc(width/2,length+2);  
  crackL = lattice*(length+2);
  lmpinput->SetAtomInfo(atominfo); //One unit cell has 4 atoms. 
  if(myrank==0)lmpinput->GenDat("eam.dat"); 
  on_crack=true;
  
  
}

void gsLMP::SetThermo(int _thermostep, double _temperature, double _pressure,string _s){
  thermoStep=_thermostep;
  temperature=_temperature;
  //rescaleStep =_scale_step;
  pressure=_pressure;
  thermoout=_s;
}

void gsLMP::SetStMin(double _percent, string _s){
  erate=_percent;
  on_erate=true;
  stressout=_s;
  on_stmin=true;
}

void gsLMP::SetStrain(double _erate, int _dstep , int _scanStep, string _s){
  erate=_erate;
  deformStep=_dstep;
  scanStep=_scanStep;
  on_erate=true;
  stressout=_s;
  on_strain=true;
}

/*void gsLMP::SetDeform(double _erate, int _dstep,int _scanStep, string _s){
  erate=_erate;
  deformStep=_dstep;
  scanStep=_scanStep;
  on_erate=true;
  stressout=_s;
  on_deform=true;
  }*/

void gsLMP::SetTSLaw(double _rcut,int _numcell,int _tsnum, string _s){
  on_tslaw=true;
  rcut=_rcut;
  numcell=_numcell;
  if(numcell==0) {
    cout <<"number of cell should be more than 2" <<endl;
    exit(1);
  }
  int size =get_global_size();
  if(numcell%size!=0) {
    cout <<"number of cell should be x # of processors" <<endl;
    exit(1);
  }

  tslawStep=_tsnum;
  tslaw=_s;

  lammps->input->one("processors 1 * 1");

}


void gsLMP::Deform(string cs){
  SetNPT(cs);
  
  delx/=timestep;

  char com[200]={NULL};
  sprintf(com,"all deform %d %s delta %f %f units box remap x",deformStep,cs.c_str(),-delx,delx);
  cmd->Fix(com);
}

/*void gsLMP::RunDeform(double estrain,string cs){

  if(thermoStep==0){ cout <<"GS: Thermo step should be set before simulation run!"<<endl; exit(1);}

  int toStep = (int)(estrain/(erate/1e12*0.001));
  cmd->Boundary("p p p");
  cmd->TimeNeighbor(timestep,2.0);
  SetComputes();

  mstep = min(thermoStep,scanStep); if(scanStep==0) mstep=thermoStep;
  int tstep = toStep/mstep+1;

  if(on_deform){
    Deform(cs);
  }else{
    char com2[200]={NULL};
    sprintf(com2,"all npt temp %f %f iso %f %f 1 drag 1",temperature,temperature,pressure,pressure);
  }
  
  cmd->Thermos(mstep,"custom step temp pe press pxx pyy pzz c_tstr c_tvol");//for output
  if(on_vmd) LMPVMD();
  if(on_cfg) LMPCFG();

  if(myrank==0) ofstress.open(stressout.c_str());
  if(myrank==0) ofthermo.open(thermoout.c_str());
  
  cmd->Run(0);
  cmd->ResetTime();
  etotal0=0;
  for(int step=0;step<tstep;step++){//total iteration! 
    PrintStress(step);
    PrintThermo(step);
    cmd->Run(mstep);
  }

  if(myrank==0)ofstress.close();
  if(myrank==0)ofthermo.close();
  }*/

void gsLMP::RunStMin(double estrain,string cs){

  if(thermoStep==0){ cout <<"GS: Thermo step should be set before simulation run!"<<endl; exit(1);}

  int toStep = (int)(estrain/erate);
  if(on_tslaw){ 
    cmd->Boundary("p s p");
  }
  else cmd->Boundary("p p p");

  cmd->TimeNeighbor(timestep,2.0);
  SetComputes();
  int tstep = toStep+1;

  if(on_stmin){
    SetNPT(cs);
  }else{
    char com2[200]={NULL};
    sprintf(com2,"all npt temp %f %f iso %f %f 1 drag 1",temperature,temperature,pressure,pressure);
  }


  mstep = min(thermoStep,scanStep); if(scanStep==0) mstep=thermoStep;
  
  cmd->Thermos(mstep,"custom step temp pe press pxx pyy pzz c_tstr c_tvol");//for output

  if(on_vmd) LMPVMD();
  if(on_cfg) LMPCFG();
  if(on_tslaw) TSLAW();

  if(myrank==0) ofstress.open(stressout.c_str());
  if(myrank==0) ofthermo.open(thermoout.c_str());
  if(myrank==0 && on_tslaw) oftslaw.open(tslaw.c_str());
  
  cmd->Run(0);
  cmd->ResetTime();
  
  delx= 0.5*lx0*erate;

  etotal0=0;
  for(int step=0;step<tstep;step++){//total iteration! 
    PrintStress(step);
    PrintThermo(step);
    

    
    char com[200]={NULL};
    sprintf(com,"change_box all %s delta %f %f remap units box",cs.c_str(),-delx,delx);  
    lammps->input->one(com);
    //cmd->MinimizeAniso(cs);
    cmd->Minimize();
    
    if(on_tslaw){
      cells->UpdateLocalRepAtom();
      cells->UpdateLocalRepStress();
    }

  }

  if(myrank==0) ofstress.close();
  if(myrank==0) ofthermo.close();
  if(myrank==0 && on_tslaw) oftslaw.close();
}

void gsLMP::RunStrain(double estrain,string cs){

  if(thermoStep==0){ cout <<"GS: Thermo step should be set before simulation run!"<<endl; exit(1);}

  int toStep = (int)(estrain/(erate/1e12*0.001));

  if(on_tslaw){
    cmd->Boundary("p s p");
  }
  else cmd->Boundary("p p p");

  cmd->TimeNeighbor(timestep,2.0);
  SetComputes();

  mstep = min(thermoStep,scanStep); if(scanStep==0) mstep=thermoStep;
  mstep = min(mstep,deformStep);

  int tstep = int(toStep/mstep)+1;

  if(on_strain){
    SetNPT(cs);
  }else{
    char com2[200]={NULL};
    sprintf(com2,"all npt temp %f %f iso %f %f 1 drag 1",temperature,temperature,pressure,pressure);
  }
  
  cmd->Thermos(mstep,"custom step temp pe press lx ly lz pxx pyy pzz c_tstr c_tvol");//for output
  if(on_vmd) LMPVMD();
  if(on_cfg) LMPCFG();
  if(on_tslaw) TSLAW();

  if(myrank==0) ofstress.open(stressout.c_str());
  if(myrank==0) ofthermo.open(thermoout.c_str());
  if(myrank==0 && on_tslaw) oftslaw.open(tslaw.c_str());
  
  cmd->Run(0);
  cmd->ResetTime();
  etotal0=0;
  
  if(on_tslaw){
    if(mstep<thermoStep) cmd->Thermos(mstep,"custom step temp pe lx ly lz press pxx pyy pzz c_tstr c_tvol");//for output
  }

  int dstep=int(deformStep/mstep);
  int tsstep=int(tslawStep/mstep);

  if(dstep==0) dstep=1;
  for(int step=0;step<tstep;step++){//total iteration! 
    PrintStress(step);
    PrintThermo(step);

    cmd->Run(mstep);

    if(on_tslaw){
      cells->UpdateLocalRepAtom();
      cells->UpdateLocalRepStress();
    }

    char com[200]={NULL};
    sprintf(com,"change_box all %s delta %f %f remap units box",cs.c_str(),-delx,delx);  
    if(step%dstep==0)lammps->input->one(com);
    if(step%20==0){
      char com[200]={NULL};
      sprintf(com,"write_restart step_%d.restart",step);
      lammps->input->one(com);
    }
  }


  if(myrank==0) ofstress.close();
  if(myrank==0) ofthermo.close();
}


void gsLMP::SetComputes(){
  cmd->Compute("csym all centro/atom fcc");
  cmd->Compute("voro all voronoi/atom");
  cmd->Compute("peatom all pe/atom");
  cmd->Compute("astres all stress/atom virial");
  cmd->Compute("tstr all reduce sum c_astres[1]");
  cmd->Compute("tvol all reduce sum c_voro[1]");
}

void gsLMP::TSLAW(){
  cells->InitCell(2,numcell,1,lammps);//cx,cy,cz..lammps
  cells->SetRcutCrackL(rcut);
  cells->UpdateLocalRepAtom();
}

void gsLMP::LMPVMD(){
    char vmdcom[200]={NULL};
    char elreplace[200]={NULL};
    sprintf(elreplace,"element %s",metal.c_str());  
    sprintf(vmdcom,"all xyz %d %s",vmdStep,vmdout.c_str());
    cmd->Dump(vmdcom);
    //cmd->DumpModi(elreplace);
}

void gsLMP::LMPCFG(){
    int natoms = lammps_get_natoms(lammps);
    char avol[200]={NULL};
    char elreplace[200]={NULL};

    sprintf(elreplace,"element %s",metal.c_str());  
    sprintf(avol,"equal vol/%d",natoms);
    cmd->Variable("avol",avol);
    
    for(int i=0;i<6;i++){
      char vname[200]={NULL};char astress[200]={NULL};
      sprintf(vname,"sig%d",i+1); sprintf(astress,"atom c_astres[%d]/10000/v_avol",i+1);
      cmd->Variable(vname,astress); //cmd->Variable("stres1","atom c_astres[1]/(20000*c_voro[1])");
    }
    
    char cfgcom[200]={NULL};
    sprintf(cfgcom,"all cfg %d %s id type xs ys zs c_peatom c_csym v_sig1 v_sig2 v_sig3",cfgStep,cfgout.c_str());
    cmd->Dump(cfgcom);
    //cmd->DumpModi(elreplace);  
}

void gsLMP::PrintStress(int step){
  myrank=get_global_rank();
  thermo=lammps->output->thermo;

  int lstep = step*mstep;

  thermo->compute_lx();
  double lx = thermo->dvalue;
  thermo->compute_ly();
  double ly = thermo->dvalue;
  thermo->compute_lz();
  double lz = thermo->dvalue;

  if(step==0){
    lx0=lx;
    ly0=ly;
    lz0=lz;
  }

  double ex = lx/lx0;
  double ey = ly/ly0;
  double ez = lz/lz0;

  thermo->compute_pxx();
  double pxx = thermo->dvalue;
  thermo->compute_pyy();
  double pyy = thermo->dvalue;
  thermo->compute_pzz();
  double pzz = thermo->dvalue; 
  if(step==0) ofstress << "step1 "<< "lx2 " << "ly3 " << "lz4 " << "ex5 " << "ey6 "<< "ez7 " <<"pxx(GPa)8 " << "pyy(GPa)9 " << "pzz(Gpa)10 "<< endl;
  if(myrank==0)ofstress << lstep<<" " <<lx <<" " <<ly<<" " << lz<<" " << ex-1.0 <<" " <<ey-1.0 <<" " <<ez-1.0 << " "<< -pxx/10000 <<" " << -pyy/10000<<" " << -pzz/10000<<endl;
}


void gsLMP::PrintThermo(int step){
  myrank=get_global_rank();
  thermo=lammps->output->thermo;

  int lstep = step*mstep;
  double dt = timestep*lstep;
  thermo->compute_temp();
  double temp = thermo->dvalue;
  thermo->compute_ke();
  double ke = thermo->dvalue;
  thermo->compute_pe();
  double pe = thermo->dvalue;
  thermo->compute_vol();
  double vol = thermo->dvalue;
  thermo->compute_press();
  double press = thermo->dvalue/10000;
  thermo->compute_enthalpy();
  double enthalpy = thermo->dvalue;
  thermo->compute_etotal();
  double etotal = thermo->dvalue;
  if(step==0)etotal0 = thermo->dvalue;
  double cons = fabs((etotal-etotal0)/etotal0);
  thermo->compute_density();
  double density = thermo->dvalue;
  
  //bars[100]kPa = 1/10000[GPa]
  if(step==0) ThermoUnit();
  if(myrank==0)ofthermo << lstep<< " " << dt << " " << temp << " " <<ke<< " "<< pe << " "<< vol <<" "<<press<<" "<< enthalpy<<" "<<etotal<<" " <<cons<< " "<<density<< endl;

}
