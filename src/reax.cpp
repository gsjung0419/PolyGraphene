//2013.11.04 G.S. Jung
#include "reax.h"
using namespace LAMMPS_NS;
gsREAX::gsREAX(){
  atoms=NULL;
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
  timestep = 0.2;
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

gsREAX::~gsREAX(){
  if(atoms!=NULL) delete atoms;
 if(lammps!=NULL)delete lammps;
  thermo=NULL;
  if(fcc!=NULL)delete fcc;
  if(cells!=NULL)delete cells;
  if(cmd!=NULL)delete cmd;
  atominfo=NULL;
  lmpinput=NULL;
}

void gsREAX::GenGraphene(double _bondl,double _sheetz,int _cx, int _cy, int _cz){
  bondl = _bondl;
  sheetz =_sheetz;
  cx = _cx;
  cy = _cy;
  cz = _cz;

  graphene->SetBondSheet(bondl, sheetz);
  atominfo = graphene->GenByCell(cx,cy,cz); //dimension..300/400/16
  atominfo->GenVMD("grahpene.xyz","C");
  lmpinput->SetAtomInfo(atominfo); //One unit cell has 4 atoms. 
  if(myrank==0)lmpinput->GenDat("reax.dat"); 

}
void gsREAX::ReadDat(){
  lammps->input->one("processors * * 1");
  char tmpeam[100]={NULL};
  sprintf(tmpeam,"%s 1 2 3 4 5 6",potential.c_str());
  cmd->ReadDataREAX("reax.dat",tmpeam);
}

void gsREAX::Init(int argc, char *argv[]){
  lmpinput=new LMPInput;
  //fcc=new FCC;
  graphene= new Graphene; 
  cells=new Cells;
  lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
  cmd=new LMPCommand;   
  thermo=lammps->output->thermo;
  cmd->REAXInit(lammps);
  myrank=get_global_rank();
  lmpinput->poType=1;//SET POTENTISL TYPE REAX
}

void gsREAX::SetAtoms(string s){
  potential=s;
  ntype=6;
  atoms = new string[6];
  atoms[0]="C";
  atoms[1]="H";
  atoms[2]="O";
  atoms[3]="N";
  atoms[4]="S";
  atoms[5]="Si";
  
}

void gsREAX::ThermoUnit(){
  if(myrank==0) ofthermo<< "step " << "time(fs) "<< "3temp(K) "<< "4ke(Kcal/mol) "<<"5pe(Kcal/mol) " << "6vol(Angstroms^3) "<< "7press(GPa) "<<"8enthalpy(Kcal/mol) " <<"9etotal(Kcal/mol) " <<"con(Et/Et0) " "11density(gram/cm^3) " <<endl;
}

void gsREAX::SetNPT(string cs){
    lx0 = lammps->domain->boxhi[0]-lammps->domain->boxlo[0];
    ly0 = lammps->domain->boxhi[1]-lammps->domain->boxlo[1];
    lz0 = lammps->domain->boxhi[2]-lammps->domain->boxlo[2];

    char com2[200]={NULL};
    string tmpx = "x";
    string tmpy = "y";
    string tmpz = "z";
    if(tmpx == cs){
      sprintf(com2,"all npt temp %f %f 1 y %f %f 1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure,pressure,pressure);
      delx = 0.5*erate/1.0e15*lx0*deformStep*timestep;
      if(on_tslaw){
	sprintf(com2,"all npt temp %f %f 1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure);
      }
    }
    else if(tmpy == cs){
      sprintf(com2,"all npt temp %f %f 1 x %f %f  1 z %f %f 1 drag 1",temperature,temperature,pressure,pressure,pressure,pressure);
      delx = 0.5*erate/1.0e15*ly0*deformStep*timestep;
    }
    else if(tmpz == cs){
      sprintf(com2,"all npt temp %f %f 1 x %f %f  1 y %f %f 1 drag 1",temperature,temperature,pressure,pressure,pressure,pressure);
      delx = 0.5*erate/1.0e15*lz0*deformStep*timestep;
    }
    cmd->Fix(com2);
  
}


void gsREAX::ReadRestart(string _readfile){
  char com[200]={NULL};
  sprintf(com,"read_restart %s",_readfile.c_str());
  
  lammps->input->one(com);
  char command2[CHNUM]={NULL};
  
  char tmpeam[100]={NULL};
  sprintf(tmpeam,"%s 1 2 3 4 5 6",potential.c_str());
  sprintf(command2,"pair_coeff  \t * * %s",tmpeam);
  lammps->input->one("pair_style  \t reax");
  lammps->input->one(command2);
}
