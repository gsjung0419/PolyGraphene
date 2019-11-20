//2013.09.28 G.S. Jung
#include "atominfo.h"

AtomInfo::AtomInfo()
{
  numAtom=0;
  box=NULL;
  subox=NULL;

  atype=NULL;
  ids=NULL;
  charge=NULL;
  x=NULL;
  y=NULL;
  z=NULL;
  vx=NULL;
  vy=NULL;
  vz=NULL;

  blockx=1;
  blocky=1;
  blockz=1;
}

AtomInfo::~AtomInfo()
{
  delete[] box;
  delete[] atype;
  delete[] ids;
  delete[] charge;
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
}

void AtomInfo::GenVMD(string _s,string _sa){
  char file[50]={NULL};
  ofstream fp;
  int ssize = _s.size();

  memcpy(file,_s.c_str(),ssize);

  fp.open(file);
  
  fp << numAtom<<endl;
  fp << "Atoms. Timestep: 0\n";
  for(int i=0;i<numAtom;i++){
    //fp<<"Al " << x[i] << " " << y[i] <<" " << z[i]<<endl;
    fp << _sa.c_str()<<" ";
    fp<< x[i] << " " << y[i] <<" " << z[i]<<endl;
  }
  fp.close();
}

void AtomInfo::SetupArray(){

  if(numAtom==0){
    cout <<"The number of atom should be setting \n";
    exit(1);
  }else{
    atype=new int[numAtom];
    ids=new int[numAtom];
    charge=new double[numAtom];
    x=new double[numAtom];
    y=new double[numAtom];
    z=new double[numAtom];

    vx=new double[numAtom];
    vy=new double[numAtom];
    vz=new double[numAtom];
    for(int i=0;i<numAtom;i++){
      atype[i]=0;
      ids[i]=0;
      x[i]=0.0;
      y[i]=0.0;
      z[i]=0.0;
      vx[i]=0.0;
      vy[i]=0.0;
      vz[i]=0.0;
    }
  }
}

void AtomInfo::SetAtomNum(int _numAtom){
  numAtom = _numAtom;
  SetupArray();
  //This automatically creates the array for the information.
}

void AtomInfo::SetBox(double * _h){
  box=new double[6]; 
  for(int i=0;i<6;i++) box[i] = _h[i];
}

void AtomInfo::SetCoord(double *_x, double *_y,double* _z)
{
  if(numAtom == 0){
    cout << "Error Set numAtom before Set Coordinates\n";
    exit(1);
  }else{
    for(int i=0;i<numAtom;i++){
      x[i]=_x[i];     
      y[i]=_y[i];     
      z[i]=_z[i];     
    }
  }
}

void AtomInfo::SetAtype(int *_atype){
  if(numAtom == 0){
    cout << "Error Set numAtom before set Atom Type\n";
    exit(1);
  }else{
    for(int i=0;i<numAtom;i++){
      atype[i] =1;//Graphene
    }
  }
}

void AtomInfo::SetAtomIds(int *_ids){
  if(numAtom == 0){
    cout << "Error Set numAtom before set Atom Type\n";
    exit(1);
  }else{
    for(int i=0;i<numAtom;i++){
      ids[i] =_ids[i];//Graphene
    }
  }
}


void AtomInfo::AssignCell(std::vector<int> &clist, AtomInfo *gatominfo){
  for(int i=0;i<numAtom;i++){
    int id = clist[i];
    x[i]=gatominfo->GetX(id);
    y[i]=gatominfo->GetY(id);
    z[i]=gatominfo->GetZ(id);
    atype[i]=gatominfo->GetAtype(id);
  }

  
}
