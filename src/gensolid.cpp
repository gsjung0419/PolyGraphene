//2013.10.01 G.S. Jung @LAMM
#include "gensolid.h"


void GenSolid::RemoveTriangle(double width, double length){
  double *h = atominfo->box;
  double midx = 0.5*(h[1]-h[0]);
  
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
    double clx = lm*y[i]+ lx-x[i];
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

void GenSolid::RemoveYTriangle(double width, double length){
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

void GenSolid::RemoveBox(double *_h)
{
  double rbox[6];
  
  for(int i=0;i<6;i++){
    rbox[i] = _h[i];
  }

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
    if(rbox[0]<=x[i] && rbox[1]>=x[i]){
      if(rbox[2]<=y[i] && rbox[3]>=y[i]){
	if(rbox[4]<=z[i] && rbox[5]>=z[i]){
	  rlist.push_back(i);
	}
      }
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

void GenSolid::MakeSpace(double dx, double dy, double dz)
{
  double *box = atominfo->GetBox();

  box[0]-=dx;
  box[1]+=dx;
  box[2]-=dy;
  box[3]+=dy;
  box[4]-=dz;
  box[5]+=dz;

}

