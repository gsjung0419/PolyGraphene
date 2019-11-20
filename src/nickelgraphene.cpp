//2013.10.03 G.S. Jung @LAMM
#include "nickelgraphene.h"

Graphene::Graphene(){
  atominfo=NULL;
  lattice=0.0;
  bondl=1.42;
  mbondl=0.0;
  sheetz=3.35;
}

Graphene::~Graphene(){
  delete atominfo;
}


AtomInfo* Graphene::GenByCell(int _lat_x, int _lat_y, int _lat_z){
  //Graphne bond length 1.42 A 
  //Sheet length 3.35 A
  int replix=_lat_x;
  int repliy=_lat_y;
  int repliz=_lat_z;


  double box[6]={0.0};

  double sqrt3 = sqrt(3);
  double latticex = sqrt3*bondl;
  double latticey = 3*bondl;
  double latticez = sheetz;
  
  double x1 = 0.0;
  double y1 = 0.5*bondl;

  double x2 = 0.5*sqrt3*bondl;
  double y2 = bondl;

  double x3 = 0.5*sqrt3*bondl;
  double y3 = 2*bondl;

  double x4 = 0.0;
  double y4 = 2.5*bondl;

  int atomNum = replix*repliy*repliz*4;
  double *x=new double[atomNum];
  double *y=new double[atomNum];
  double *z=new double[atomNum];
  int *atype=new int[atomNum];
  int *ids=new int[atomNum];
  
  for(int i=0;i<replix;i++){
    for(int j=0;j<repliy;j++){
      for(int k=0;k<repliz;k++){
	double x0 = i*latticex;
	double y0 = j*latticey;
	double z0 = (k+0.5)*latticez;
	
	int index = 4*(i*repliy*repliz + j*repliz + k);

	x[index] = x0+x1;
	y[index] = y0+y1;
	z[index] = z0;

	x[index+1] = x0+x2;
	y[index+1] = y0+y2;
	z[index+1] = z0;

	x[index+2] = x0+x3;
	y[index+2] = y0+y3;
	z[index+2] = z0;

	x[index+3] = x0+x4;
	y[index+3] = y0+y4;
	z[index+3] = z0;
      
      }
    }
  }

  for(int i=0;i<atomNum;i++){
    atype[i]=1;
    ids[i]=1+i;
    //cout << x[i] << " "<<y[i] <<" " << z[i]<<endl;
  }
  
  if(atominfo!=NULL){
    cout<< "AtomInfo is not NULL pointer\n";
  }else atominfo=new AtomInfo;

  atominfo->SetAtomNum(atomNum);
  box[0]=0;
  box[1]=replix*latticex;
  box[2]=0.0;
  box[3]=repliy*latticey;
  box[4]=0.0;
  box[5]=repliz*latticez;

  atominfo->SetBox(box);
  atominfo->SetCoord(x,y,z);
  atominfo->SetAtype(atype);
  atominfo->SetAtomIds(ids);
  
  return atominfo;
}

AtomInfo* Graphene::GenByDimension(double lx, double ly, double lz){
  //Graphne bond length 1.42 A 
  //Sheet length 3.35 A

  double sqrt3 = sqrt(3);
  double latticex = sqrt3*bondl;
  double latticey = 3*bondl;
  double latticez = sheetz;  
  
  int replix= int(lx/latticex);
  int repliy= int(ly/latticey);
  int repliz= int(lz/latticez);
  
  if(replix%2 != 0) replix++;

  return GenByCell(replix,repliy,1);

}

void Graphene::GenGB(double theta){
  //Graphne bond length 1.42 A 
  //Sheet length 3.35 A
  double *box = atominfo->box;

  double *x = atominfo->x;
  double *y = atominfo->y;
  double *z = atominfo->z;
  
  double midx = 0.5*(box[1]);
  int natom = atominfo->numAtom;
  for(int i=0;i<natom;i++){
    if(x[i]<midx){

      double nx,ny;
      double cx = x[i]-midx;

      nx = cx*cos(theta)+y[i]*sin(theta);
      ny = -cx*sin(theta)+y[i]*cos(theta);

      double siftxx = -0.5*box[1]*cos(theta);
      double siftxy = 0.5*box[1]*sin(theta);
      
      double siftyx = -box[3]*sin(theta);
      double siftyy = -box[3]*cos(theta);

      nx+=midx;

      if(nx>=midx){ nx+=siftxx; ny+=siftxy;}
      if(ny>box[3]) {nx+=siftyx; ny+=siftyy;}


      x[i]=nx;
      y[i]=ny;

    }else if(x[i]>=midx){
      double nx,ny;
      double cx = x[i]-midx;
      nx = cx*cos(theta)-y[i]*sin(theta);
      ny = cx*sin(theta)+y[i]*cos(theta);

      double siftxx = 0.5*box[1]*cos(theta);
      double siftxy = 0.5*box[1]*sin(theta);
      
      double siftyx = box[3]*sin(theta);
      double siftyy = -box[3]*cos(theta);

      nx+=midx;

      if(nx<midx){ nx+=siftxx; ny+=siftxy;}
      if(ny>box[3]) {nx+=siftyx; ny+=siftyy;}

      x[i]=nx;
      y[i]=ny;

    }
  }

  double rbox[6];

  rbox[0] = box[1];
  rbox[1] = 2*box[1];
  rbox[2] = box[2];
  rbox[3] = box[3];
  rbox[4] = box[4];
  rbox[5] = box[5];

  RemoveBox(rbox);

  rbox[0] = -box[1];
  rbox[1] = 0;
  rbox[2] = box[2];
  rbox[3] = box[3];
  rbox[4] = box[4];
  rbox[5] = box[5];

  RemoveBox(rbox);

}


void Graphene::RemoveCenterx(double ycrackl, double width)
{
  double rbox[6];
  double *box = atominfo->GetBox();

  double xlength = box[1]-box[0];
  double sqrt3 = sqrt(3);
  //double width = sqrt3*num_wA*bondl;

  rbox[0] = box[0]+0.5*xlength-0.5*width;
  rbox[1] = box[0]+0.5*xlength+0.5*width;
  rbox[2] = box[2];
  rbox[3] = box[2]+ycrackl;
  rbox[4] = box[4];
  rbox[5] = box[5];

  RemoveBox(rbox);

}
