//2013.11.04 G.S. Jung
#include "eam.h"
using namespace LAMMPS_NS;

gsEAM::gsEAM(){
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
  timestep = 0.001;
  deformStep=0;
  thermoStep=0;
  rescaleStep=0;
  outputStep=0;
  scanStep=0;
  crackL=0.0;
  lattice=0.0;
  erate=0.0;
  lx0=0.0;
  ly0=0.0;
  lz0=0.0;

  on_crack=false;
  on_deform=false;
  on_strain=false;
  on_stmin=false;
  on_vmd=false;
  on_cfg=false;
  }

gsEAM::~gsEAM(){
  if(lammps!=NULL)delete lammps;
  lammps=NULL;
  thermo=NULL;
  if(fcc!=NULL)delete fcc;
  if(cells!=NULL)delete cells;
  if(cmd!=NULL)delete cmd;
  atominfo=NULL;
  lmpinput=NULL;
  }

void gsEAM::GenFCC(double _lc,double _rcut,int _cx, int _cy, int _cz){
  lattice = _lc;
  rcut =_rcut;
  cx = _cx;
  cy = _cy;
  cz = _cz;

  fcc->SetLattice(lattice);
  atominfo = fcc->GenByCell(cx,cy,cz); //dimension..300/400/16
  atominfo->GenVMD("fcc.xyz","Al");


  lmpinput->SetAtomInfo(atominfo); //One unit cell has 4 atoms. 
  if(myrank==0)lmpinput->GenDat("eam.dat"); 

}

void gsEAM::ReadDat(){
  lammps->input->one("processors 1 * 1");
  char tmpeam[100]={NULL};
  sprintf(tmpeam,"%s %s",potential.c_str(),metal.c_str());
  cmd->ReadDataEAM("eam.dat",tmpeam);
}

void gsEAM::GenFCC110(double _lc,double _rcut,int _cx, int _cy, int _cz){
  lattice = _lc;
  rcut =_rcut;
  cx = _cx;
  cy = _cy;
  cz = _cz;

  fcc->SetLattice(lattice);
  atominfo = fcc->GenByCell110(cx,cy,cz); //dimension..300/400/16
  atominfo->GenVMD("fcc110.xyz","Al");

  lmpinput->SetAtomInfo(atominfo); //One unit cell has 4 atoms. 
  if(myrank==0)lmpinput->GenDat("eam.dat"); 
}

void gsEAM::Init(int argc, char *argv[]){
  lmpinput=new LMPInput;
  fcc=new FCC;
  cells=new Cells;
  lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
  cmd=new LMPCommand;   
  thermo=lammps->output->thermo;
  cmd->EAMInit(lammps);
  myrank=get_global_rank();

}

void gsEAM::SetMetal(string s1,string s2){
  if(metal==""){ metal=s1;
  }
  else{cout <<"The metal name is already set!!"<<endl; exit(1); }
  potential=s2;
}

void gsEAM::Run(int toStep,int ostep){
  outputStep=ostep;
  if(on_crack){
    if(!on_erate){ cout <<"Strain rate should be set before simulation run!"<<endl; exit(1);}
    if(deformStep==0){ cout <<"Deform step should be set before simulation run!"<<endl; exit(1);}
    if(scanStep==0){ cout <<"Scan step should be set before simulation run!"<<endl; exit(1);}
  }
  
  if(thermoStep==0){ cout <<"Thermo step should be set before simulation run!"<<endl; exit(1);}
  if(outputStep==0){ cout <<"Output step should be set before simulation run!"<<endl; exit(1);}

  lammps->input->one("processors 1 * 1");
  char tmpeam[100]={NULL};
  sprintf(tmpeam,"%s %s",potential.c_str(),metal.c_str());
  //cmd->ReadDataEAM("eam.dat",potential.c_str());
  cmd->ReadDataEAM("eam.dat",tmpeam);
  cmd->MinimizeISO();
  lammps->input->one("change_box all boundary p f p");
  cmd->Minimize();
  cmd->TimeNeighbor(timestep,2.0);
  cmd->ResetTime();

  cmd->Compute("csym all centro/atom fcc");
  cmd->Compute("voro all voronoi/atom");
  cmd->Compute("peatom all pe/atom");
  cmd->Compute("astres all stress/atom virial");
  cmd->Compute("tstr all reduce sum c_astres[1]");
  cmd->Compute("tvol all reduce sum c_voro[1]");
  //cmd->Variable("tstres", "equal tst");
  //cmd->Variable("tvol", "equal tvol");
  //cmd->Variable("tmp","equal v_stres/v_tvol");

  cells->InitCell(2,numcell,1,lammps);
  cells->SetRcutCrackL(rcut);
  cells->UpdateLocalRepAtom();
  
  double lx = lammps->domain->boxhi[0]-lammps->domain->boxlo[0];
  delx = 0.5*erate/1.0e12*lx*deformStep*timestep;
  
  char command0[200]={NULL};
  sprintf(command0,"change_box all x delta %f %f remap units box",-delx,delx);  
  lammps->input->one(command0);


  cmd->Thermos(scanStep,"custom step lx ly lz temp press c_tstr c_tvol");
  cmd->Fix("all nve");

  char command[200]={NULL};
  if(rescaleStep!=0){
    sprintf(command,"all temp/rescale %d 0.0 0.0 0.0 1.0",rescaleStep);
    cmd->Fix(command);
  }else{
    sprintf(command,"all temp/rescale %d 0.0 0.0 0.0 1.0",rescaleStep);
    //cmd->Fix("all temp/rescale 1 0.0 0.0 0.0 1.0");
  }

  
  for(int i=0;i<6;i++){
    char vname[200]={NULL};char astress[200]={NULL};
    sprintf(vname,"sig%d",i+1); sprintf(astress,"atom c_astres[%d]/20000",i+1);
    cmd->Variable(vname,astress); //cmd->Variable("stres1","atom c_astres[1]/(20000*c_voro[1])");
  }
  

  char command1[200]={NULL};
  sprintf(command1,"all xyz %d %s_Result.xyz",outputStep,metal.c_str());
  cmd->Dump(command1);

  char command2[200]={NULL};
  sprintf(command2,"element %s",metal.c_str());
  cmd->DumpModi(command2);

  if(on_cfg){
    char command3[200]={NULL};
    sprintf(command3,"all cfg %d %s*.cfg id type xs ys zs c_peatom ",cfgStep,cfgout.c_str());
    cmd->Dump(command3);
    cmd->DumpModi(command2);
  }


  cmd->Run(scanStep);

  int dstep=deformStep/scanStep;
  for(int step=1;step<toStep/scanStep;step++){//total iteration! 
    if(step%dstep==0){
      //lammps->input->one("change_box all x delta -1.0 1.0 remap units box");
      char command[200]={NULL};
      sprintf(command,"change_box all x delta %f %f remap units box",-delx,delx);  
      lammps->input->one(command0);
    }

    cells->UpdateLocalRepAtom();
    cells->UpdateLocalRepStress();
    cmd->Run(scanStep);
  }
}

void gsEAM::ThermoUnit(){
  if(myrank==0) ofthermo<< "step " << "time(ps) "<< "3temp(K) "<< "4ke(eV) "<<"5pe(eV) " << "6vol(Angstroms^3) "<< "7press(GPa) "<<"8enthalpy(eV) " <<"9etotal(eV) " <<"con(Et/Et0) " "11density(gram/cm^3) " <<endl;
  
}

void gsEAM::SetNPT(string cs){
    lx0 = lammps->domain->boxhi[0]-lammps->domain->boxlo[0];
    ly0 = lammps->domain->boxhi[1]-lammps->domain->boxlo[1];
    lz0 = lammps->domain->boxhi[2]-lammps->domain->boxlo[2];
    double lyh =lammps->domain->boxhi[1];
    //double lyl =lyh-4.1;
    double lyl =lammps->domain->boxlo[1];
    double lxh =lammps->domain->boxhi[0];
    double lxl =lammps->domain->boxlo[0];
    double lzh =lammps->domain->boxhi[2];
    double lzl =lammps->domain->boxlo[2];

    char com2[200]={NULL};
    string tmpx = "x";
    string tmpy = "y";
    string tmpz = "z";
    

    if(tmpx == cs){
      sprintf(com2,"all npt temp %f %f 1 y %f %f 1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure,pressure,pressure);
      //cmd->Fix("all nve");      
      delx = 0.5*erate/1.0e12*lx0*deformStep*timestep;
      if(on_tslaw){
	char com[200]={NULL};
	//sprintf(com,"region yedge block INF INF %f %f INF INF",lyl,lyh);
	sprintf(com,"region yedge block INF INF INF INF INF INF");
	lammps->input->one(com);
	lammps->input->one("group ywall region yedge");
	lammps->input->one("fix contem ywall temp/rescale 1 0 0 0 1.0");

	sprintf(com2,"all npt temp %f %f 1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure);
      }
    }
    else if(tmpy == cs){
      sprintf(com2,"all npt temp %f %f 1 x %f %f  1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure,pressure,pressure);
      delx = 0.5*erate/1.0e12*ly0*deformStep*timestep;
    }
    else if(tmpz == cs){
      sprintf(com2,"all npt temp %f %f 1 x %f %f  1 y %f %f 1 drag 1",temperature,temperature,pressure,pressure,pressure,pressure);
      delx = 0.5*erate/1.0e12*lz0*deformStep*timestep;
    }
    cmd->Fix(com2);
  
}

void gsEAM::ReadRestart(string _readfile){
  char com[200]={NULL};
  sprintf(com,"read_restart %s",_readfile.c_str());
  
  lammps->input->one(com);
  char command2[CHNUM]={NULL};
  
  char tmpeam[100]={NULL};
  sprintf(tmpeam,"%s %s",potential.c_str(),metal.c_str());
  sprintf(command2,"pair_coeff  \t * * %s",tmpeam);
  lammps->input->one("pair_style  \t eam/alloy");
  lammps->input->one(command2);
}
