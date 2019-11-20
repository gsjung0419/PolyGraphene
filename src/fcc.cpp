//2013.10.03 G.S. Jung @LAMM
#include "fcc.h"

FCC::FCC(){
  atominfo=NULL;
  lattice=0.0;
  bondl=0.0;
  mbondl=0.0;
  sheetz=0.0;
}

FCC::~FCC(){
  delete atominfo;
}


void FCC::SetLattice(double _lattice)
{
  lattice = _lattice;
}

AtomInfo* FCC::GenByCell(int _lat_x, int _lat_y, int _lat_z){
  int replix=_lat_x;
  int repliy=_lat_y;
  int repliz=_lat_z;

  double base = 0.25*lattice;

  double box[6]={0.0};

  double x1 = base;
  double y1 = base;
  double z1 = base;

  double x2 = base+0.5*lattice;
  double y2 = base+0.5*lattice;
  double z2 = base+0.0;

  double x3 = base+0.5*lattice;
  double y3 = base;
  double z3 = base+0.5*lattice;

  double x4 = base;
  double y4 = base+0.5*lattice;
  double z4 = base+0.5*lattice;

  
  int atomNum = replix*repliy*repliz*4;

  double *x=new double[atomNum];
  double *y=new double[atomNum];
  double *z=new double[atomNum];
  int *atype=new int[atomNum];
  int *ids=new int[atomNum];
  
  for(int i=0;i<replix;i++){
    for(int j=0;j<repliy;j++){
      for(int k=0;k<repliz;k++){
	double x0 = i*lattice;
	double y0 = j*lattice;
	double z0 = k*lattice;
	
	int index = 4*(i*repliy*repliz + j*repliz + k);

	x[index] = x0+x1;
	y[index] = y0+y1;
	z[index] = z0+z1;

	x[index+1] = x0+x2;
	y[index+1] = y0+y2;
	z[index+1] = z0+z2;

	x[index+2] = x0+x3;
	y[index+2] = y0+y3;
	z[index+2] = z0+z3;

	x[index+3] = x0+x4;
	y[index+3] = y0+y4;
	z[index+3] = z0+z4;
      
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
  box[1]=replix*lattice;
  box[2]=0.0;
  box[3]=repliy*lattice;
  box[4]=0.0;
  box[5]=repliz*lattice;

  atominfo->SetBox(box);
  atominfo->SetCoord(x,y,z);
  atominfo->SetAtype(atype);
  atominfo->SetAtomIds(ids);
  
  return atominfo;
}

AtomInfo* FCC::GenByCell110(int _lat_x, int _lat_y, int _lat_z){
  int replix=_lat_x;
  int repliy=_lat_y;
  int repliz=_lat_z;

  double root2 = sqrt(2.0);
  double lattice_x = lattice/root2;
  double lattice_y = lattice/root2;
  double lattice_z = lattice;

  double base_x = 0.25*lattice_x;
  double base_y = 0.25*lattice_y;
  double base_z = 0.25*lattice_z;

  double box[6]={0.0};

  double x1 = base_x;
  double y1 = base_y;
  double z1 = base_z;

  double x2 = base_x+0.5*lattice_x;
  double y2 = base_y+0.5*lattice_y;
  double z2 = base_z+0.5*lattice_z;

  int atomNum = replix*repliy*repliz*2;

  double *x=new double[atomNum];
  double *y=new double[atomNum];
  double *z=new double[atomNum];
  int *atype=new int[atomNum];
  int *ids=new int[atomNum];
  
  for(int i=0;i<replix;i++){
    for(int j=0;j<repliy;j++){
      for(int k=0;k<repliz;k++){
	double x0 = i*lattice_x;
	double y0 = j*lattice_y;
	double z0 = k*lattice_z;
	
	int index = 2*(i*repliy*repliz + j*repliz + k);

	x[index] = x0+x1;
	y[index] = y0+y1;
	z[index] = z0+z1;

	x[index+1] = x0+x2;
	y[index+1] = y0+y2;
	z[index+1] = z0+z2;

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
  box[1]=replix*lattice_x;
  box[2]=0.0;
  box[3]=repliy*lattice_y;
  box[4]=0.0;
  box[5]=repliz*lattice_z;

  atominfo->SetBox(box);
  atominfo->SetCoord(x,y,z);
  atominfo->SetAtype(atype);
  atominfo->SetAtomIds(ids);

  cout <<"Atom number of the system " << atomNum<<endl;
  
  return atominfo;
}

AtomInfo* FCC::GenByCell111(int _lat_x, int _lat_y, int _lat_z){
  int replix=_lat_x;
  int repliy=_lat_y;
  int repliz=_lat_z;

  double root2 = sqrt(2.0);
  double lattice_x = lattice/root2;
  double lattice_y = lattice/root2;
  double lattice_z = lattice;

  double base_x = 0.25*lattice_x;
  double base_y = 0.25*lattice_y;
  double base_z = 0.25*lattice_z;

  double box[6]={0.0};

  double x1 = base_x;
  double y1 = base_y;
  double z1 = base_z;

  double x2 = base_x+0.5*lattice_x;
  double y2 = base_y+0.5*lattice_y;
  double z2 = base_z+0.5*lattice_z;

  int atomNum = replix*repliy*repliz*2;

  double *x=new double[atomNum];
  double *y=new double[atomNum];
  double *z=new double[atomNum];
  int *atype=new int[atomNum];
  int *ids=new int[atomNum];
  
  for(int i=0;i<replix;i++){
    for(int j=0;j<repliy;j++){
      for(int k=0;k<repliz;k++){
	double x0 = i*lattice_x;
	double y0 = j*lattice_y;
	double z0 = k*lattice_z;
	
	int index = 2*(i*repliy*repliz + j*repliz + k);

	x[index] = x0+x1;
	y[index] = y0+y1;
	z[index] = z0+z1;

	x[index+1] = x0+x2;
	y[index+1] = y0+y2;
	z[index+1] = z0+z2;

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
  box[1]=replix*lattice;
  box[2]=0.0;
  box[3]=repliy*lattice;
  box[4]=0.0;
  box[5]=repliz*lattice;

  atominfo->SetBox(box);
  atominfo->SetCoord(x,y,z);
  atominfo->SetAtype(atype);
  atominfo->SetAtomIds(ids);

  cout <<"Atom number of the system " << atomNum<<endl;
  
  return atominfo;
}

AtomInfo* FCC::GenByDimension(double lx, double ly, double lz){
  
  int replix= int(lx/lattice);
  int repliy= int(ly/lattice);
  int repliz= int(lz/lattice);

  AtomInfo *atominfo = GenByCell(replix,repliy,repliz);

  return atominfo;

}

void FCC::RemoveCenterx(double ycrackl, double width)
{
  double rbox[6];
  double *box = atominfo->GetBox();

  double xlength = box[1]-box[0];


  rbox[0] = box[0]+0.5*xlength-0.5*width;
  rbox[1] = box[0]+0.5*xlength+0.5*width;
  rbox[2] = box[2];
  rbox[3] = box[2]+ycrackl;
  rbox[4] = box[4];
  rbox[5] = box[5];

  RemoveBox(rbox);
  //for(int i=0;i<6;i++)cout <<rbox[i]<<endl;

}

void FCC::RemoveCenterxc(int x_width,int y_length)
{
  //cout <<"start removecenterxc"<<endl;
  double rbox[6];
  double *box = atominfo->GetBox();
  
  double xlength = box[1]-box[0];
  int total_x = int((box[1]-box[0])/lattice);
  int x_position = total_x/2;
  //cout <<"x_position " <<x_position<<endl;

  if(total_x%2 !=0){ 
    printf("The cell is not x-symmetric!!\n");
    exit(1);
  }

  int width = x_width/2;
  //cout <<"width " <<x_position<<endl;


  if(x_width!=1){
    if(x_width%2!=0){
      printf("The center removed cell is not symmetric!!\n");
      exit(1);
    }
    rbox[0] = box[0]+lattice*(x_position-width);
    rbox[1] = box[0]+lattice*(x_position+width); 
    rbox[2] = box[2];
    rbox[3] = box[2]+lattice*y_length;
    rbox[4] = box[4];
    rbox[5] = box[5];
  }else{
    rbox[0] = box[0]+lattice*((double)x_position-0.5*x_width);
    rbox[1] = box[0]+lattice*((double)x_position+0.5*x_width); 
    rbox[2] = box[2];
    rbox[3] = box[2]+lattice*y_length;
    rbox[4] = box[4];
    rbox[5] = box[5];

  }


  //cout <<"end removecenterxc"<<endl;
  RemoveBox(rbox);

  //for(int i=0;i<6;i++)cout <<rbox[i]<<endl;

}
