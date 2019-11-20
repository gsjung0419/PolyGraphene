//2013.09.27 G.S. Jung LAMM
#include "readfile.h"


void ReadFile::ReadCoord(string s, int tstep)
{
  ifstream fp;
  int timestep=tstep;
  char filename[100]={NULL};
  memcpy(filename,s.c_str(),s.size());

  fp.open(filename);
  atominfo=new AtomInfo;
  char dumchar[100];

  int numAtom=0;
  int ftimestep=0;

  fp>>numAtom;
  atominfo->SetAtomNum(numAtom);//Set the Number of Atom
  string tmp1,tmp2, tmp3;
  fp>>tmp1>>tmp2>>tmp3;
  cout <<tmp1 << " " <<tmp2<<" "<< tmp3<<endl;
  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();
  int *atype =atominfo->GetAtype();
  int *ids = atominfo->ids;

  double xmax=0;
  double xmin=0;
  double ymax=0;
  double ymin=0;
  double zmax=0;
  double zmin=0;
  double dx=2.0;
  double dy=2.0;
  double dz=2.0;


  for(int i=0;i<numAtom;i++){

    //getline(fp,line);

    string tmp;
    //istringstream(line)>>atype[i]>>x[i]>>y[i]>>z[i];
    //istringstream(line)>>x[i]>>y[i]>>z[i];
    fp>>tmp>>x[i]>>y[i]>>z[i];


    ids[i] = i+1;
    if(i== 0){
      xmax=x[i];
      xmin=x[i];
      ymax=y[i];
      ymin=y[i];
      zmax=z[i];
      zmin=z[i];
    }

    double tdx=fabs(x[i]-x[0]);
    double tdy=fabs(y[i]-y[0]);
    double tdz=fabs(z[i]-z[0]);

    if(xmax < x[i]) xmax=x[i];
    if(xmin > x[i]) xmin=x[i];
    if(ymax < y[i]) ymax=y[i];
    if(ymin > y[i]) ymin=y[i];
    if(zmax < z[i]) zmax=z[i];
    if(zmin > z[i]) zmin=z[i];


    //For debugging coordinate
    atype[i]=1;
    //cout << ids[i]<<" " << atype[i] << " 0 " <<x[i] <<" " <<y[i] <<" " <<z[i]<<endl;    

  }

  double box[6];
  box[0] = xmin-0.5*dx;
  box[1] = xmax+0.5*dx;
  box[2] = ymin-0.5*dy;
  box[3] = ymax+0.5*dy;
  box[4] = zmin-0.5*dz;
  box[5] = zmax+0.5*dz;

  atominfo->SetBox(box);

  double xlo=box[0];
  double xhi=box[1];
  double ylo=box[2];
  double yhi=box[3];
  double zlo=box[4];
  double zhi=box[5];

  double lengthy = (yhi-ylo)*0.5;
  ylo -=lengthy;
  yhi +=lengthy;


  /*cout <<"LAMMPS DATA FILE\n\n"; 
  cout <<numAtom <<" atoms\n\n";

  cout <<"2 atom types\n\n";

  cout <<xlo<<" " <<xhi<<" " <<"xlo xhi\n";
  cout <<ylo<<" "  <<yhi<<" " <<"ylo yhi\n";
  cout <<"-67 67 zlo zhi\n\n";
  cout <<"Atoms";*/

  fp.close();
}

void ReadFile::ReadPDB(string s)
{
  ifstream fp;
  char filename[100]={NULL};
  memcpy(filename,s.c_str(),s.size());

  fp.open(filename);
  atominfo=new AtomInfo;
  char dumchar[100];

  int numAtom=0;
  int ftimestep=0;
  string line;

  fp >> numAtom;
  atominfo->SetAtomNum(numAtom);//Set the Number of Atom

  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();
  int *atype =atominfo->GetAtype();
  int *ids = atominfo->ids;

  double xmax=0;
  double xmin=0;
  double ymax=0;
  double ymin=0;
  double zmax=0;
  double zmin=0;
  double dx=2.0;
  double dy=2.0;
  double dz=2.0;

  for(int i=0;i<numAtom;i++){

    string tmp1, atom, res, num, tmp2,tmp3,tmp4;
    int id;
    double tx,ty,tz;
    fp>>tmp1>>id>>atom>>res>>num>>tx>>ty>>tz>>tmp2>>tmp3>>tmp4;
    ids[i]=id;
    x[i]=tx;
    y[i]=ty;
    z[i]=tz;

    //ids[i] = i+1;
    if(i== 0){
      xmax=x[i];
      xmin=x[i];
      ymax=y[i];
      ymin=y[i];
      zmax=z[i];
      zmin=z[i];
    }

    double tdx=fabs(x[i]-x[0]);
    double tdy=fabs(y[i]-y[0]);
    double tdz=fabs(z[i]-z[0]);

    if(xmax < x[i]) xmax=x[i];
    if(xmin > x[i]) xmin=x[i];
    if(ymax < y[i]) ymax=y[i];
    if(ymin > y[i]) ymin=y[i];
    if(zmax < z[i]) zmax=z[i];
    if(zmin > z[i]) zmin=z[i];


    //For debugging coordinate
    atype[i]=1;
    //cout << ids[i]<<" " << atype[i] << " 0 " <<x[i] <<" " <<y[i] <<" " <<z[i]<<endl;    

  }

  double box[6];
  box[0] = xmin-0.5*dx;
  box[1] = xmax+0.5*dx;
  box[2] = ymin-0.5*dy;
  box[3] = ymax+0.5*dy;
  box[4] = zmin-0.5*dz;
  box[5] = zmax+0.5*dz;

  atominfo->SetBox(box);

  double xlo=box[0];
  double xhi=box[1];
  double ylo=box[2];
  double yhi=box[3];
  double zlo=box[4];
  double zhi=box[5];
  
  double lengthy = (yhi-ylo)*0.5;
  ylo -=lengthy;
  yhi +=lengthy;


  cout <<"LAMMPS DATA FILE\n\n"; 
  cout <<numAtom <<" atoms\n\n";

  cout <<"2 atom types\n\n";

  cout <<xlo<<" " <<xhi<<" " <<"xlo xhi\n";
  cout <<ylo<<" "  <<yhi<<" " <<"ylo yhi\n";
  cout <<"-67 67 zlo zhi\n\n";
  cout <<"Atoms";

  fp.close();
}


void  ReadFile::ReadBox(string s){

  ifstream fp;
  char filename[100]={NULL};
  memcpy(filename,s.c_str(),s.size());  
  
  fp.open(filename);
  double box[6];
  fp >> box[0] >> box[1] >> box[2] >>box[3] >> box[4] >> box[5];
  atominfo->SetBox(box);

  fp.close();
}

void ReadFile::ReadCrack(string s, int nstep)
{
  ifstream fp;
  char filename[100]={NULL};
  memcpy(filename,s.c_str(),s.size());
  cout << s <<endl;
  fp.open(filename,ios::in);
  char dumchar[100];

  int numAtom=0;
  string line;

  /*vector <double> x1;
  vector <double> y1;
  vector <double> z1;
  vector <double> x2;
  vector <double> y2;
  vector <double> z2;*/

  double *x1 = new double[nstep];
  double *y1 = new double[nstep];
  double *z1 = new double[nstep];
  double *stres11 = new double[nstep];
  double *stres12 = new double[nstep];
  double *stres13 = new double[nstep];
  double *stres14 = new double[nstep];
  double *stres15 = new double[nstep];
  double *stres16 = new double[nstep];

  double *x2 = new double[nstep];
  double *y2 = new double[nstep];
  double *z2 = new double[nstep];
  double *stres21 = new double[nstep];
  double *stres22 = new double[nstep];
  double *stres23 = new double[nstep];
  double *stres24 = new double[nstep];
  double *stres25 = new double[nstep];
  double *stres26 = new double[nstep];

  

  for(int i=0;i<nstep;i++){
    x1[i]=0.0;
    y1[i]=0.0;
    z1[i]=0.0;
    x2[i]=0.0;
    y2[i]=0.0;
    z2[i]=0.0;
  }

  int aid1;
  int aid2;  
  for(int istep=0;istep<nstep;istep++){
    for(int i=0;i<9;i++){ 
      getline(fp,line);
      //cout<<line<<endl;
    }
      getline(fp,line);
    istringstream(line) >> aid1 >> x1[istep] >>y1[istep] >>z1[istep]>>stres11[istep]>>stres12[istep]>>stres13[istep];
    getline(fp,line);
    istringstream(line) >> aid2 >> x2[istep] >>y2[istep] >>z2[istep]>>stres21[istep]>>stres22[istep]>>stres23[istep];
  }


  for(int istep=0;istep<nstep;istep++){
    cout << 0.1*fabs(x1[istep] - x2[istep])<< " " << 0.5*(stres11[istep]+stres21[istep])<<endl;
  }

  fp.close();
}

void ReadFile::RemoveTriangle(double width, double length){
  double *h = atominfo->box;
  double midx = 0.5*(h[1]+h[0]);
  
  double hm = -0.5*width/length;
  double lm = 0.5*width/length;

  double hx = midx+0.5*width;
  double lx = midx-0.5*width;

  vector <int>rlist;
  
  int numAtom = atominfo->GetNumAtom();
  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();
  double *vx = atominfo->GetVX();
  double *vy = atominfo->GetVY();
  double *vz = atominfo->GetVZ();
  int *atype = atominfo->GetAtype();
  int *Ids = atominfo->GetAtomIds();
  double *charge = atominfo->GetCharge();

  for(int i=0;i<numAtom;i++){
    double chx = hm*y[i]+ hx -x[i];
    double clx = lm*y[i]+ lx -x[i];
    if(chx>0 && clx<0){
      rlist.push_back(i);
    }
  }  
  
  int size_rlist = rlist.size();
  //check the number of removed atom! 
  //cout << "size_rlist " << size_rlist <<" Atom Num "<< numAtom <<endl;

  int count=0;
  for(int i=0;i<size_rlist;i++){
    int rindex=rlist[i]-count;
    int ncopy = numAtom-rindex-1;
    memmove(&x[rindex],&x[rindex+1],ncopy*sizeof(double));
    memmove(&y[rindex],&y[rindex+1],ncopy*sizeof(double));
    memmove(&z[rindex],&z[rindex+1],ncopy*sizeof(double));
    memmove(&atype[rindex],&atype[rindex+1],ncopy*sizeof(int));
    memmove(&Ids[rindex],&Ids[rindex+1],ncopy*sizeof(int));
    memmove(&charge[rindex],&charge[rindex+1],ncopy*sizeof(double));
    memmove(&vx[rindex],&vx[rindex+1],ncopy*sizeof(double));
    memmove(&vy[rindex],&vy[rindex+1],ncopy*sizeof(double));
    memmove(&vz[rindex],&vz[rindex+1],ncopy*sizeof(double));

    numAtom--;
    count++;
  }
  
  atominfo->ReAtomNum(numAtom);

}

void ReadFile::RemoveYTriangle(double width, double length){
  double *h = atominfo->box;
  double midy = 0.5*(h[3]-h[2]);
  
  double hm = -0.5*width/length;
  double lm = 0.5*width/length;

  double hx = midy+0.5*width;
  double lx = midy-0.5*width;

  vector <int>rlist;
  
  int numAtom = atominfo->GetNumAtom();
  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();
  double *vx = atominfo->GetVX();
  double *vy = atominfo->GetVY();
  double *vz = atominfo->GetVZ();
  int *atype = atominfo->GetAtype();
  int *Ids = atominfo->GetAtomIds();
  double *charge = atominfo->GetCharge();

  for(int i=0;i<numAtom;i++){
    double chy = hm*x[i]+ hx -y[i];
    double cly = lm*x[i]+ lx-y[i];
    if(chy>0 && cly<0){
      rlist.push_back(i);
    }
  }  
  
  int size_rlist = rlist.size();
  //check the number of removed atom! 
  //cout << "size_rlist " << size_rlist <<" Atom Num "<< numAtom <<endl;

  int count=0;
  for(int i=0;i<size_rlist;i++){
    int rindex=rlist[i]-count;
    int ncopy = numAtom-rindex-1;
    memmove(&x[rindex],&x[rindex+1],ncopy*sizeof(double));
    memmove(&y[rindex],&y[rindex+1],ncopy*sizeof(double));
    memmove(&z[rindex],&z[rindex+1],ncopy*sizeof(double));
    memmove(&atype[rindex],&atype[rindex+1],ncopy*sizeof(int));
    memmove(&Ids[rindex],&Ids[rindex+1],ncopy*sizeof(int));
    memmove(&charge[rindex],&charge[rindex+1],ncopy*sizeof(double));
    memmove(&vx[rindex],&vx[rindex+1],ncopy*sizeof(double));
    memmove(&vy[rindex],&vy[rindex+1],ncopy*sizeof(double));
    memmove(&vz[rindex],&vz[rindex+1],ncopy*sizeof(double));

    numAtom--;
    count++;
  }
  
  atominfo->ReAtomNum(numAtom);

}


void ReadFile::Translate(double dx, double dy,double dz){
  int numAtom = atominfo->GetNumAtom();
  double *x = atominfo->GetX();
  double *y = atominfo->GetY();
  double *z = atominfo->GetZ();

  for(int i=0;i<numAtom;i++){
    x[i]+=dx;
    y[i]+=dy;
    z[i]+=dz;
    
  }

}
