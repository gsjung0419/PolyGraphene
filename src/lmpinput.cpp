//2013.10.07 G.S. Jung @LAMM
#include "lmpinput.h"

LMPInput::LMPInput(){
  atominfo=NULL;
  poType=0;
  numAtomType=1;
  snum=1;
}

LMPInput::~LMPInput(){
  atominfo=NULL;
}

void LMPInput::GenDatTypeNum(string s, int tnum){
 char filename[100]={NULL};
  memcpy(filename,s.c_str(),s.size());
  ofstream fp;
  fp.open(filename);

  int numAtom = atominfo->numAtom;
  int numType = atominfo->numType;
  int *atype = atominfo->atype;
  int *ids = atominfo->ids;
  double *x = atominfo->x;
  double *y = atominfo->y;
  double *z = atominfo->z;
  double *box =atominfo->box;

  double xlo,xhi;
  double ylo,yhi;
  double zlo,zhi;

  xlo=box[0];
  xhi=box[1];
  ylo=box[2];
  yhi=box[3];
  zlo=box[4];
  zhi=box[5];

  fp<<"LAMMPS DATA FILE\n\n"; 
  fp<<numAtom <<" atoms\n\n";


  if(poType==EAM)fp<<numAtomType<<" atom types\n\n";
  if(poType==REAX)fp<<"6"<<" atom types\n\n";

  
  fp<<xlo<<" " <<xhi<<" " <<"xlo xhi\n";
  fp<<ylo<<" " <<yhi<<" " <<"ylo yhi\n";
  fp<<zlo<<" " <<zhi<<" " <<"zlo zhi\n\n";

  /*if(poType==REAX){
    fp<<"Masses\n\n";
    fp<<"1 12.01\n";
    fp<<"2 1.001\n";
    fp<<"3 15.9994\n";
    fp<<"4 14.0067\n";
    fp<<"5 32.065\n";
    fp<<"6 28.0855\n\n";
    }*/

  fp<<"Atoms\n\n";

  for(int i=0;i<numAtom;i++){
    if(poType==EAM) fp<<i+snum<<" "<< tnum <<" " << x[i] << " " << y[i] <<" " << z[i]<<endl;
    //if(poType==REAX) fp<<ids[i]<<" "<< atype[i] << " 0  " << x[i] << " " << y[i] <<" " << z[i]<<endl;
    if(poType==REAX) fp<< i+snum<<" "<< tnum << " 0  " << x[i] << " " << y[i] <<" " << z[i]<<endl;

  }
  fp.close();
}

void LMPInput::GenDat(string s){
  char filename[100]={NULL};
  memcpy(filename,s.c_str(),s.size());
  ofstream fp;
  fp.open(filename);

  int numAtom = atominfo->numAtom;
  int numType = atominfo->numType;
  int *atype = atominfo->atype;
  int *ids = atominfo->ids;
  double *x = atominfo->x;
  double *y = atominfo->y;
  double *z = atominfo->z;
  double *box =atominfo->box;

  double xlo,xhi;
  double ylo,yhi;
  double zlo,zhi;

  xlo=box[0];
  xhi=box[1];
  ylo=box[2];
  yhi=box[3];
  zlo=box[4];
  zhi=box[5];

  fp<<"LAMMPS DATA FILE\n\n"; 
  fp<<numAtom <<" atoms\n\n";


  if(poType==EAM)fp<<numAtomType<<" atom types\n\n";
  if(poType==REAX)fp<<"2"<<" atom types\n\n";

  
  fp<<xlo<<" " <<xhi<<" " <<"xlo xhi\n";
  fp<<ylo<<" " <<yhi<<" " <<"ylo yhi\n";
  fp<<zlo<<" " <<zhi<<" " <<"zlo zhi\n\n";

  /*if(poType==REAX){
    fp<<"Masses\n\n";
    fp<<"1 12.01\n";
    fp<<"2 1.001\n";
    fp<<"3 15.9994\n";
    fp<<"4 14.0067\n";
    fp<<"5 32.065\n";
    fp<<"6 28.0855\n\n";
    }*/

  fp<<"Atoms\n\n";

  for(int i=0;i<numAtom;i++){
    if(poType==EAM) fp<<i+snum<<" "<< atype[i] <<" " << x[i] << " " << y[i] <<" " << z[i]<<endl;
    //if(poType==REAX) fp<<ids[i]<<" "<< atype[i] << " 0  " << x[i] << " " << y[i] <<" " << z[i]<<endl;
    if(poType==REAX) fp<< i+snum<<" "<< atype[i] << " 0  " << x[i] << " " << y[i] <<" " << z[i]<<endl;

  }
  fp.close();
  
}

void LMPInput::GenInput(string s){
  char filename[50]={NULL};
  memcpy(filename,s.c_str(),s.size());
  ofstream fp;
  fp.open(filename);
  
  int numAtom = atominfo->GetNumAtom();
  int numType = atominfo->GetNumType();
  int *atype = atominfo->GetAtype();
  int *ids = atominfo->GetAtomIds();
  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();
  double *box =atominfo->GetBox();

  double xlo,xhi;
  double ylo,yhi;
  double zlo,zhi;

  xlo=box[0];
  xhi=box[1];
  ylo=box[2];
  yhi=box[3];
  zlo=box[4];
  zhi=box[5];
    

  fp<<"LAMMPS DATA FILE\n\n"; 
  fp<<numAtom <<" atoms\n\n";

  fp<<"6 atom types\n\n";
  
  fp<<xlo<<" " <<xhi<<" " <<"xlo xhi\n";
  fp<<ylo<<" " <<yhi<<" " <<"ylo yhi\n";
  fp<<zlo<<" " <<zhi<<" " <<"zlo zhi\n\n";

  fp<<"Masses\n\n";

  fp<<"1 12.01\n";
  fp<<"2 1.001\n";
  fp<<"3 15.9994\n";
  fp<<"4 14.0067\n";
  fp<<"5 32.065\n";
  fp<<"6 28.0855\n\n";

  fp<<"Atoms\n\n";

  for(int i=0;i<numAtom;i++){
    fp<<ids[i]<<" "<< atype[i] << " 0  " << x[i] << " " << y[i] <<" " << z[i]<<endl;
    //fp<<i+1<<" "<< atype[i] <<" " << x[i] << " " << y[i] <<" " << z[i]<<endl;
  }
  fp.close();

}

void LMPInput::GenVMD(string _s){
  char file[50]={NULL};
  ofstream fp;
  int ssize = _s.size();

  memcpy(file,_s.c_str(),ssize);

  fp.open(file);
  
  int numAtom = atominfo->GetNumAtom();
  int *atype = atominfo->GetAtype();
  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();
  double *box =atominfo->GetBox();

  double xlo,xhi;
  double ylo,yhi;
  double zlo,zhi;

  xlo=box[0];
  xhi=box[1];
  ylo=box[2];
  yhi=box[3];
  zlo=box[4];
  zhi=box[5];
  
  fp << numAtom<<endl;
  fp << "Atoms. Timestep: 0\n";
  for(int i=0;i<numAtom;i++){
    fp<<"C " << x[i] << " " << y[i] <<" " << z[i]<<endl;
  }
  fp.close();


}



