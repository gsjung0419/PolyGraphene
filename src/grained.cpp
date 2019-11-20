//2014.06.06 G.S. Jung @LAMM
#include "grained.h"

Grained::Grained(){
  atominfo=NULL;
  blist=NULL;
  clist=NULL;
  chead=NULL;
  sav_x=NULL;
  sav_y=NULL;
  rotate=NULL;
  particles = 72;
  x_min = -180, x_max = 180;
  y_min = -360, y_max = 360;
  z_min = -50, z_max = 50;
  c_x=60,c_y=120;
  maxhanum=112000;

  /*particles = 4;
  x_min = -50, x_max = 50;
  y_min = -50, y_max = 50;
  z_min = -50, z_max = 50;
  c_x=20,c_y=20;
  maxhanum=4000;

  particles = 20;
  x_min = -100, x_max = 100;
  y_min = -100, y_max = 100;
  z_min = -50, z_max = 50;
  c_x=40,c_y=40;
  maxhanum=16000;*/

  n_x=25,n_y=25,n_z=5;

  ibondcut2=4.0;
  addcut2=1.2*1.2;

  myrank = get_global_rank();
  size = get_global_size();

  hfname="history.xyz";
  char filename[100]={NULL};
  memcpy(filename,hfname.c_str(),hfname.size());
  hfp.open(filename);
  hstep=0;
 
  mover_x = x_max;
  srand(time(NULL));

  //Cell List Initialize
  

}

Grained::~Grained(){
  if(atominfo!=NULL)delete atominfo;
}

void Grained::GenVoro(){
  Graphene *gra;
  AtomInfo *atominfo;

  gra = new Graphene;
  gra->SetBondSheet(1.42,3.35);

  using namespace voro;
  const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

  //atominfo = gra->GenByDimension(2.0*(x_max-x_min),2.0*(y_max-y_min),1);
  atominfo = gra->GenByDimension(2.0*(y_max-y_min),2.0*(y_max-y_min),1);
  int numAtom = atominfo->numAtom;

  double *x,*y,*z;

  x = atominfo->x;
  y = atominfo->y;
  z = atominfo->z;

  for(int i=0;i<numAtom;i++){
    //x[i] -= x_max*2.0;
    x[i] -= y_max*2.0;
    y[i] -= y_max*2.0;
  }

  if(myrank==0){
    atominfo->GenVMD("graphene.xyz","C");
  }
  
  //container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
  con = new container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);

  ifstream fip;
  char filename[100]={NULL};
  string s="points_p.gnu";
  memcpy(filename,s.c_str(),s.size());
  
  fip.open(filename);
  
  /*for(int i=0;i<particles;i++) {
    int tmpindx;
    double tmpx,tmpy,tmpz;
    tmpx=x_min+rnd()*(x_max-x_min);
    tmpy=y_min+rnd()*(y_max-y_min);
    tmpz=0.0;
    con->put(i,tmpx,tmpy,tmpz);
    } */ //random generation of voronoi seeds 

  //int nnx=18;
  //int nny=36;
  
  for(int i=0;i<particles;i++) {
    int tmpindx;
    double tmpx,tmpy,tmpz;
    tmpx=x_min+rnd()*(x_max-x_min);
    tmpy=y_min+rnd()*(y_max-y_min);
    tmpz=0.0;
    con->put(i,tmpx,tmpy,tmpz);
  } //random generation of voronoi seeds 
  
  if(rotate==NULL) rotate = new double[particles];
  for(int i=0;i<particles;i++){
    rotate[i]=0.0;
  }

  for(int i=0;i<particles;i++){
    rotate[i]=rnd();
    if(myrank==0)cout <<i <<" rotate " << rotate[i] <<endl;
  }

  global_dsum(rotate,particles);
  atomnum=0;

  //if(myrank==0){
  for(int ci = 0; ci<numAtom;ci++){
    int id;
    double rx,ry,rz;
    double r_x,r_y,r_z;
    for(int pi =0;pi<particles;pi++){
      double theta = 3.141592*rotate[pi]*0.5;
      double rsin = sin(theta);
      double rcos = cos(theta);
      r_x = x[ci]*rcos-y[ci]*rsin;
      r_y = x[ci]*rsin+y[ci]*rcos;
      //if(con.find_voronoi_cell(r_x,r_y,z[ci],rx,ry,rz,id))
      if(con->find_voronoi_cell(r_x,r_y,z[ci],rx,ry,rz,id))
	if(id==pi){
	  atomnum++;
	  gx.push_back(r_x);
	  gy.push_back(r_y);
	  gz.push_back(r_z);
	}
    }
  }
  //}

  BuildBlist();//build bond list and checklist

  if(myrank==0){
    FILE *f1=safe_fopen("grained.xyz","w");
    fprintf(f1,"%d\n",atomnum);
    fprintf(f1,"Atoms. Timestep: 0\n");
    for(int i=0;i<gx.size();i++){
      if(check[i]) fprintf(f1,"C %f %f %f\n",gx[i],gy[i],0.0);    
    }

    fclose(f1);

    SaveCoord();
    
    FILE *f3=safe_fopen("grained0.xyz","w");
    fprintf(f3,"%d\n",gx.size());
    fprintf(f3,"Atoms. Timestep: 0\n");
    for(int i=0;i<gx.size();i++){
      fprintf(f3,"C %f %f %f\n",gx[i],gy[i],gz[i]);    
    }
    
    fclose(f3);

  }

  GenLMP();
  
  delete gra;

  gx.clear();
  gy.clear();
  gz.clear();

  if(myrank==0){
    double vvol=con->sum_cell_volumes();
    printf("Container volume : %g\n"
	   "Voronoi volume   : %g\n"
	   "Difference       : %g\n",cvol,vvol,vvol-cvol);
    
    // Output the particle positions in gnuplot format
    con->draw_particles("points_p.gnu");
    
    // Output the Voronoi cells in gnuplot format
    con->draw_cells_gnuplot("points_v.gnu");  
  }
  
  
}


void Grained::GenVoroRe(){
  Graphene *gra;
  AtomInfo *atominfo;

  gra = new Graphene;
  gra->SetBondSheet(1.42,3.35);

  using namespace voro;
  const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

  //atominfo = gra->GenByDimension(2.0*(x_max-x_min),2.0*(y_max-y_min),1);
  atominfo = gra->GenByDimension(2.0*(y_max-y_min),2.0*(y_max-y_min),1);
  int numAtom = atominfo->numAtom;

  double *x,*y,*z;

  x = atominfo->x;
  y = atominfo->y;
  z = atominfo->z;

  for(int i=0;i<numAtom;i++){
    //x[i] -= x_max*2.0;
    x[i] -= y_max*2.0;
    y[i] -= y_max*2.0;
  }

  if(myrank==0){
    atominfo->GenVMD("graphene.xyz","C");
  }
  
  //container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
  con = new container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);

  ifstream fip;
  char filename[100]={NULL};
  string s="points_p.gnu";
  memcpy(filename,s.c_str(),s.size());
  
  fip.open(filename);


  /*for(int i=0;i<particles;i++) {
    int tmpindx;
    double tmpx,tmpy,tmpz;
    tmpx=x_min+rnd()*(x_max-x_min);
    tmpy=y_min+rnd()*(y_max-y_min);
    
    tmpz=0.0;
    con->put(i,tmpx,tmpy,tmpz);
    }*///random generation of voronoi seeds 

  for(int i=0;i<particles;i++) {
    int tmpindx;
    double tmpx,tmpy,tmpz;
    fip>>tmpindx >> tmpx >> tmpy >> tmpz;
    
    tmpz=0.0;
    con->put(i,tmpx,tmpy,tmpz);
  }
  
  if(rotate==NULL) rotate = new double[particles];
  for(int i=0;i<particles;i++){
    rotate[i]=0.0;
  }

  for(int i=0;i<particles;i++){
    rotate[i]=rnd();
    if(myrank==0)cout <<i <<" rotate " << rotate[i] <<endl;
  }

  global_dsum(rotate,particles);
  atomnum=0;

  //if(myrank==0){
  for(int ci = 0; ci<numAtom;ci++){
    int id;
    double rx,ry,rz;
    double r_x,r_y,r_z;
    for(int pi =0;pi<particles;pi++){
      double theta = 3.141592*rotate[pi]/3.0;
      double rsin = sin(theta);
      double rcos = cos(theta);
      r_x = x[ci]*rcos-y[ci]*rsin;
      r_y = x[ci]*rsin+y[ci]*rcos;
      //if(con.find_voronoi_cell(r_x,r_y,z[ci],rx,ry,rz,id))
      if(con->find_voronoi_cell(r_x,r_y,z[ci],rx,ry,rz,id))
	if(id==pi){
	  atomnum++;
	  gx.push_back(r_x);
	  gy.push_back(r_y);
	  gz.push_back(r_z);
	}
    }
  }
  //}

  //BuildBlist();//build bond list and checklist

  if(myrank==0){
    FILE *f1=safe_fopen("grained.xyz","w");
    fprintf(f1,"%d\n",atomnum);
    fprintf(f1,"Atoms. Timestep: 0\n");
    for(int i=0;i<gx.size();i++){
      //if(check[i]) fprintf(f1,"C %f %f %f\n",gx[i],gy[i],0.0);    
    }

    fclose(f1);

    //SaveCoord();
    
    FILE *f3=safe_fopen("grained0.xyz","w");
    fprintf(f3,"%d\n",gx.size());
    fprintf(f3,"Atoms. Timestep: 0\n");
    for(int i=0;i<gx.size();i++){
      //fprintf(f3,"C %f %f %f\n",gx[i],gy[i],gz[i]);    
    }
    
    fclose(f3);

  }

  //GenLMP();
  
  delete gra;

  gx.clear();
  gy.clear();
  gz.clear();

  if(myrank==0){
    double vvol=con->sum_cell_volumes();
    printf("Container volume : %g\n"
	   "Voronoi volume   : %g\n"
	   "Difference       : %g\n",cvol,vvol,vvol-cvol);
    
    // Output the particle positions in gnuplot format
    con->draw_particles("points_p.gnu");
    
    // Output the Voronoi cells in gnuplot format
    con->draw_cells_gnuplot("points_v.gnu");  
  }


}

void Grained::GenLMP(){
  using namespace voro;

  if(myrank==0){
    FILE *f2=safe_fopen("grained.data","w");
    fprintf(f2,"LAMMPS DATA FILE\n\n");
    fprintf(f2,"%d atoms\n\n",atomnum);
    fprintf(f2,"2 atom types\n\n");
    fprintf(f2,"%f %f xlo xhi\n",x_min-20,x_max+20);
    fprintf(f2,"%f %f ylo yhi\n",y_min-20,y_max+20);
    fprintf(f2,"%f %f zlo zhi\n\n",z_min,z_max);
    
    fprintf(f2,"Atoms\n\n");
    int id_c = 0;
    
    for(int i=0;i<gx.size();i++){
      if(check[i]){ 
	id_c++;
	fprintf(f2,"%d 1 0 %f %f %f\n",id_c,gx[i],gy[i],gz[i]*0.5);    
      }
    }
    
    fclose(f2);

    id_c=0;
    hfp<<maxhanum<<endl;
    hfp<<"Atoms. Timestep: "<<hstep<<endl;

    int *cgrain=new int[particles];
    for(int i=0;i<particles;i++) cgrain[i]=0;

    double *agrain=new double[particles];
    for(int i=0;i<particles;i++) agrain[i]=0;

    for(int i=0;i<gx.size();i++){
      if(check[i]){ 
	id_c++;
	//fprintf(f2,"%d 1 0 %f %f %f\n",id_c,gx[i],gy[i],gz[i]*0.5);    
	double sinv=0.0;
	double deg60 = 3.141592/3.0;
	for(int j =0;j<blist[i].size();j++){ 
	  int ij=blist[i][j];
	  double dx = (gx[ij]-gx[i]);
	  double dy = (gy[ij]-gy[i]);
	  double dr = sqrt(dx*dx+dy*dy);
	  if(dx>0.0) sinv = asin(dy/dr);	  
	}

	if(sinv>0.0) while(sinv>=deg60) sinv-=deg60;
	if(sinv<0.0) while(sinv<=0.0) sinv+=deg60;

	double r_x = gx[i];
	double r_y = gy[i];
	double rx,ry,rz;
	int id;
	//if(con.find_voronoi_cell(r_x,r_y,0.0,rx,ry,rz,id))
	con->find_voronoi_cell(r_x,r_y,0.0,rx,ry,rz,id);
	cgrain[id]++;
	agrain[id]+=sinv;
	//hfp << "C " << gx[i] <<" " << gy[i]<< " " << 0.01*rotate[id] <<endl;
	

	//hfp << "C " << gx[i] <<" " << gy[i]<< " " << sinv/deg60*0.01 <<endl;
	//hfp << "C " << gx[i] <<" " << gy[i]<< " " << 0.1*sinv/deg60-0.05 <<endl;
      
      }
    }

    id_c=0;

    for(int i=0;i<gx.size();i++){
      if(check[i]){ 
	id_c++;
	//fprintf(f2,"%d 1 0 %f %f %f\n",id_c,gx[i],gy[i],gz[i]*0.5);    
	double sinv=0.0;
	double deg60 = 3.141592/3.0;
	for(int j =0;j<blist[i].size();j++){ 
	  int ij=blist[i][j];
	  double dx = (gx[ij]-gx[i]);
	  double dy = (gy[ij]-gy[i]);
	  double dr = sqrt(dx*dx+dy*dy);
	  if(dx>0.0) sinv = asin(dy/dr);	  
	}

	if(sinv>0.0) while(sinv>=deg60) sinv-=deg60;
	if(sinv<0.0) while(sinv<=0.0) sinv+=deg60;

	double r_x = gx[i];
	double r_y = gy[i];
	double rx,ry,rz;
	int id;
	//if(con.find_voronoi_cell(r_x,r_y,0.0,rx,ry,rz,id))
	con->find_voronoi_cell(r_x,r_y,0.0,rx,ry,rz,id);
	//grain[id]++;
	//agrain[id]+=sinv;
	//hfp << "C " << gx[i] <<" " << gy[i]<< " " << 0.01*rotate[id] <<endl;

	//hfp << "C " << gx[i] <<" " << gy[i]<< " " << 0.01*agrain[id]/cgrain[id] <<endl;
	hfp << "C " << gx[i] <<" " << gy[i]<< " " << 0.01*rotate[id] <<endl;
      
      }
    }

    for(int i=0;i<maxhanum-id_c;i++){
      hfp<< "C " << -100+x_min << -100+x_min << " 0.0" <<endl;
    }
    
    hstep++;

    delete[] cgrain;
    delete[] agrain;

  }//End of myrank==0
  


}

void Grained::BuildBlist(){
  
  if(blist==NULL)blist = new vector<int>[gx.size()];
  else{
    for(int i=0;i<tatom;i++){
      blist[i].clear();
    }

    delete[] blist;
    blist=NULL;
  }

  if(check==NULL)check = new bool[gx.size()];
  else{
    delete[] check;
    check=NULL;
  }

  tatom =gx.size();

  BuildClist();

  if(blist==NULL)blist = new vector<int>[tatom];
  if(check==NULL)check = new bool[tatom];
  for(int i=0;i<tatom;i++){
    check[i]=false;
  }

  for(int cell_i=0;cell_i<c_x;cell_i++){
    for(int cell_j=0;cell_j<c_y;cell_j++){
      int nca = clist[cell_i][cell_j].size();
      
      for(int cai=0;cai<nca;cai++){
	int i=clist[cell_i][cell_j][cai];
	
	for(int ci=cell_i-1;ci<=cell_i+1;ci++){
	  for(int cj=cell_j-1;cj<=cell_j+1;cj++){
	    if(ci>-1 && cj>-1 && ci<c_x && cj<c_y){
	      for(int ca =0;ca<clist[ci][cj].size();ca++){
		int ja=clist[ci][cj][ca];
		if(i!=ja){
		  double dx = gx[i]-gx[ja];
		  double dy = gy[i]-gy[ja];
		  double dz = gz[i]-gz[ja];
		  
		  double dist = dx*dx + dy*dy + dz*dz;
		  if(dist < ibondcut2) {
		    blist[i].push_back(ja);
		  }	      
		}
	      }
	    }
	  }
	}
      }
    }
  }//Build Bond List using Cell list
  
  atomnum=0;
    for(int i=0;i<tatom;i++){
    if(blist[i].size()==3 && gx[i]<x_max  && gx[i]>x_min){
      atomnum++;
      check[i]=true;
    }
    if(blist[i].size()==2&& gx[i]<x_max  && gx[i]>x_min){
      int i1 = blist[i][0];
      int i2 = blist[i][1];
      if(blist[i1].size()==3 && blist[i2].size()==3){
	atomnum++;
      	check[i]=true;
      }
    }
  }

}

void Grained::BuildBlist2(){
  if(blist==NULL)blist = new vector<int>[gx.size()];
  else{
    for(int i=0;i<tatom;i++){
      blist[i].clear();
    }

    delete[] blist;
    blist=NULL;
  }

  if(check==NULL)check = new bool[gx.size()];
  else{
    delete[] check;
    check=NULL;
  }

  if(myrank==0)tatom =gx.size();
  else{tatom=0;}
  
  global_isum(&tatom,1);

  double *x,*y,*z;

  x = new double[tatom];
  y = new double[tatom];
  z = new double[tatom];
  
  if(myrank==0){
    for(int i=0;i<tatom;i++){
      x[i]=gx[i];
      y[i]=gy[i];
      z[i]=gz[i];
    }
  }else{
    for(int i=0;i<tatom;i++){
      x[i]=0.0;
      y[i]=0.0;
      z[i]=0.0;
    }
  }
  
  global_dsum(x,tatom);
  global_dsum(y,tatom);
  global_dsum(z,tatom);

  BuildClist();

  if(blist==NULL)blist = new vector<int>[tatom];
  if(check==NULL)check = new bool[tatom];
  for(int i=0;i<tatom;i++){
    check[i]=false;
  }

  for(int cell_i=0;cell_i<c_x;cell_i++){
    for(int cell_j=0;cell_j<c_y;cell_j++){
      int nca = clist[cell_i][cell_j].size();
      
      for(int cai=0;cai<nca;cai++){
	int i=clist[cell_i][cell_j][cai];
	
	for(int ci=cell_i-1;ci<=cell_i+1;ci++){
	  for(int cj=cell_j-1;cj<=cell_j+1;cj++){
	    if(ci>-1 && cj>-1 && ci<c_x && cj<c_y){
	      for(int ca =0;ca<clist[ci][cj].size();ca++){
		int ja=clist[ci][cj][ca];
		if(i!=ja){
		  double dx = gx[i]-gx[ja];
		  double dy = gy[i]-gy[ja];
		  double dz = gz[i]-gz[ja];
		  
		  double dist = dx*dx + dy*dy + dz*dz;
		  if(dist < ibondcut2) {
		    blist[i].push_back(ja);
		  }	      
		}
	      }
	    }
	  }
	}
      }
    }
  }//Build Bond List using Cell list

  atomnum=0;

  int countd=0;

  for(int i=0;i<tatom;i++){
    if(blist[i].size()==3){
      atomnum++;
      check[i]=true;
      
      int i1 = blist[i][0];
      int i2 = blist[i][1];
      int i3 = blist[i][2];
      
      double rx1 = gx[i1]-gx[i];
      double ry1 = gy[i1]-gy[i];

      double rx2 = gx[i2]-gx[i];
      double ry2 = gy[i2]-gy[i];

      double rx3 = gx[i3]-gx[i];
      double ry3 = gy[i3]-gy[i];
     
      double dot,srx1,srx2,srx3;
      dot = rx1*rx2+ry1*ry2;
      srx1 = sqrt(rx1*rx1 + ry1*ry1);
      srx2 = sqrt(rx2*rx2 + ry2*ry2);
      
      dot=dot/srx1/srx2;

      if((dot<-0.85 || dot>0.0) && check[i]){
	if(myrank==0)cout <<"degrees check!! Cosine value " << dot<< endl;
	atomnum--;
	check[i]=false;
	countd++;
      }

      dot = rx1*rx3+ry1*ry3;
      srx1 = sqrt(rx1*rx1 + ry1*ry1);
      srx3 = sqrt(rx3*rx3 + ry3*ry3);
      
      dot=dot/srx1/srx3;

      if((dot<-0.85 || dot>0.0) && check[i]){
	if(myrank==0)cout <<"60 degrees check!! Cosine value " << dot<< endl;
	atomnum--;
	check[i]=false;
	countd++;
      }

      dot = rx2*rx3+ry2*ry3;
      srx2 = sqrt(rx2*rx2 + ry2*ry2);
      srx3 = sqrt(rx3*rx3 + ry3*ry3);
      
      dot=dot/srx2/srx3;

      if((dot<-0.85 || dot>0.0) && check[i]){
	if(myrank==0)cout <<"60 degrees check!! Cosine value " << dot<< endl;
	atomnum--;
	check[i]=false;
	countd++;
      }
      
      
    }
    else if(blist[i].size()==2){
      int i1 = blist[i][0];
      int i2 = blist[i][1];
      
      if(blist[i1].size()==3 && blist[i2].size()==3){
	atomnum++;
      	check[i]=true;
      }

      double rx1 = gx[i1]-gx[i];
      double ry1 = gy[i1]-gy[i];

      double rx2 = gx[i2]-gx[i];
      double ry2 = gy[i2]-gy[i];
     
      double dot = rx1*rx2+ry1*ry2;
      double srx1 = sqrt(rx1*rx1 + ry1*ry1);
      double srx2 = sqrt(rx2*rx2 + ry2*ry2);
      
      dot=dot/srx1/srx2;

      if(dot<-0.85 && check[i]){
	if(myrank==0)cout <<"Line atom check!! Cosine value " << dot<< endl;
	atomnum--;
	check[i]=false;
      }

    }
    else if(blist[i].size()==1){
	atomnum++;
      	check[i]=true;
    }
  }

  delete []x;
  delete []y;
  delete []z;

}

void Grained::Readxyz1(){

    gx.clear();
    gy.clear();
    gz.clear();

    double *x,*y;
    if(myrank==0){    
      tatom=sav_anum;      
    }else{
      tatom=0;
    }

    global_isum(&tatom,1);

    x = new double[tatom];
    y = new double[tatom];
    
    if(myrank==0){    
      for(int i=0;i<tatom;i++){
	x[i]=sav_x[i];
	y[i]=sav_y[i];
      }
    }else{
      for(int i=0;i<tatom;i++){
	x[i]=0.0;
	y[i]=0.0;
      }
    }
    
    global_dsum(x,tatom);
    global_dsum(y,tatom);
    
    for(int i=0;i<tatom;i++){
      gx.push_back(x[i]);
      gy.push_back(y[i]);
      gz.push_back(0.0);
    }
    
    BuildBlist2();//build bond list and checklist
    
    delete[] x;
    delete[] y;

}

void Grained::Readxyz2(){
    gx.clear();
    gy.clear();
    gz.clear();
    
    double *x,*y;
    double *tx,*ty;
    ReadFile *readf;    
    readf = new ReadFile;

    if(myrank==0){    
      
      readf->ReadCoord("grained.xyz",0);
      atominfo = readf->GetAtomInfo();
    
      cout << "read size " << atominfo->numAtom<<endl;

      tx = atominfo->x;
      ty = atominfo->y;
      
      tatom = atominfo->numAtom;

    }else{tatom=0;}

    global_isum(&tatom,1);

    x = new double[tatom];
    y = new double[tatom];
    
    if(myrank==0){
      for(int i=0;i<tatom;i++){
	x[i]=tx[i];
	y[i]=ty[i];
      }
    }else{
      for(int i=0;i<tatom;i++){
	x[i]=0.0;
	y[i]=0.0;
      }
    }

    delete readf;


    global_dsum(x,tatom);
    global_dsum(y,tatom);

    for(int i=0;i<tatom;i++){
      gx.push_back(x[i]);
      gy.push_back(y[i]);
      gz.push_back(0.0);
    }


    BuildBlist3();

    if(myrank==0){    
      SaveCoord();

      /*using namespace voro;
      
      FILE *f1=safe_fopen("grained.xyz","w");
      fprintf(f1,"%d\n",gx.size());
      fprintf(f1,"Atoms. Timestep: 0\n");
      for(int i=0;i<gx.size();i++){
	fprintf(f1,"C %f %f %f\n",gx[i],gy[i],0.0);    
      }
      
      fclose(f1);*/
    }


  delete[] x;
  delete[] y;


}

void Grained::BuildBlist3(){
  if(blist==NULL)blist = new vector<int>[gx.size()];
  else{
    for(int i=0;i<tatom;i++){
      blist[i].clear();
    }
    delete[] blist;
    blist=NULL;
  }

  if(check==NULL)check = new bool[gx.size()];
  else{
    delete[] check;
    check=NULL;
  }

  tatom =gx.size();

  BuildClist();

  if(blist==NULL)blist = new vector<int>[tatom];
  if(check==NULL)check = new bool[tatom];



  for(int i=0;i<tatom;i++){
    check[i]=false;
  }

  
  if(myrank==0)cout << "Making Bond List....\n";   

  for(int cell_i=0;cell_i<c_x;cell_i++){
    for(int cell_j=0;cell_j<c_y;cell_j++){
      int nca = clist[cell_i][cell_j].size();
      
      for(int cai=0;cai<nca;cai++){
	int i=clist[cell_i][cell_j][cai];
	
	for(int ci=cell_i-1;ci<=cell_i+1;ci++){
	  for(int cj=cell_j-1;cj<=cell_j+1;cj++){
	    if(ci>-1 && cj>-1 && ci<c_x && cj<c_y){
	      for(int ca =0;ca<clist[ci][cj].size();ca++){
		int ja=clist[ci][cj][ca];
		if(i!=ja){
		  double dx = gx[i]-gx[ja];
		  double dy = gy[i]-gy[ja];
		  double dz = gz[i]-gz[ja];
		  
		  double dist = dx*dx + dy*dy + dz*dz;
		  if(dist < ibondcut2) {
		    blist[i].push_back(ja);
		  }	      
		}
	      }
	    }
	  }
	}
      }
    }
  }//Build Bond List using Cell list

  for(int i=0;i<tatom;i++){
    if(blist[i].size()==2){
      int i1 = blist[i][0];
      int i2 = blist[i][1];

      double vari = rnd();

      
      double addx = (vari+0.5)*((gx[i]-gx[i1])+(gx[i]-gx[i2]))+gx[i];
      vari =rnd();
      double addx1 = (vari+0.5)*(gx[i]-gx[i1])+gx[i];
      vari =rnd();
      double addx2 = (vari+0.5)*(gx[i]-gx[i2])+gx[i];
      vari =rnd();
      double addy = (vari+0.5)*((gy[i]-gy[i1])+(gy[i]-gy[i2]))+gy[i];
      vari =rnd();
      double addy1 = (vari+0.5)*(gy[i]-gy[i1])+gy[i];
      vari =rnd();
      double addy2 = (vari+0.5)*(gy[i]-gy[i2])+gy[i];

      //cout <<"addx " << addx << " addy "<<addy<<endl;

      bool okadd =true;
      bool okadd1 =true;
      bool okadd2 =true;
      
      for(int j=0;j<gx.size();j++){
	double dx = addx - gx[j];
	double dy = addy - gy[j];
	double cdist = dx*dx + dy*dy;

	double dx1 = addx1 - gx[j];
	double dy1 = addy1 - gy[j];
	double cdist1 = dx1*dx1 + dy1*dy1;

	double dx2 = addx2 - gx[j];
	double dy2 = addy2 - gy[j];
	double cdist2 = dx2*dx2 + dy2*dy2;

	rnd_cut();
	if(cdist <addcut2) okadd = false;
	if(cdist1 <(1.65*1.65)) okadd1 = false;
	if(cdist2 <(1.65*1.65)) okadd2 = false;
      }

      if(addx>x_max || addx<x_min || addy>y_max || addy<y_min) okadd =false;
      if(addx1>x_max || addx1<x_min || addy1>y_max || addy1<y_min) okadd1 =false;
      if(addx2>x_max || addx2<x_min || addy2>y_max || addy2<y_min) okadd2 =false;
      
      if(okadd){
	gx.push_back(addx);
	gy.push_back(addy);
	gz.push_back(0.0);
	//cout << "atom added!! \n";
      }
      else if(okadd1){
	gx.push_back(addx1);
	gy.push_back(addy1);
	gz.push_back(0.0);
	//cout << "atom added!! \n";
      }
      else if(okadd2){
	gx.push_back(addx2);
	gy.push_back(addy2);
	gz.push_back(0.0);
	//cout << "atom added!! \n";
      }


    }else if(blist[i].size()==1){
      int i1 = blist[i][0];
      double rx=(gx[i]-gx[i1]);
      double ry=(gy[i]-gy[i1]);
      double vari=rnd();
      double addx = (vari*0.2+0.9)*(0.5*rx-1.732/2.0*ry)+gx[i];
      vari=rnd();
      double addy = (vari*0.2+0.9)*(1.732/2.0*rx + 0.5*ry)+gy[i];      

      //cout <<"addx " << addx << " addy "<<addy<<endl;
      bool okadd =true;
      
      for(int j=0;j<gx.size();j++){
	double dx = addx - gx[j];
	double dy = addy - gy[j];
	double cdist = dx*dx + dy*dy;
	rnd_cut();
	if(cdist <addcut2) okadd = false;
      }
      if(addx>x_max || addx<x_min || addy>y_max || addy<y_min) okadd =false;

      //if(addx < mover_x)okadd=false;
      if(okadd){
	gx.push_back(addx);
	gy.push_back(addy);
	gz.push_back(0.0);
	//cout << "atom added!! \n";
      }

      okadd =true;
      vari=rnd();
      addx = (vari*0.2+0.9)*(0.5*rx+1.732/2.0*ry)+gx[i];
      vari=rnd();
      addy = (vari*0.2+0.9)*(-1.732/2.0*rx + 0.5*ry)+gy[i];      

      for(int j=0;j<gx.size();j++){
	double dx = addx - gx[j];
	double dy = addy - gy[j];
	double cdist = dx*dx + dy*dy;
	rnd_cut();
	if(cdist <addcut2) okadd = false;
      }

      if(addx>x_max || addx<x_min || addy>y_max || addy<y_min) okadd =false;

      //if(addx < mover_x)okadd=false;
      if(okadd){
	gx.push_back(addx);
	gy.push_back(addy);
	gz.push_back(0.0);
	//cout << "atom added!! \n";
      }
      
    }
  }

  delete[] check;
  check=NULL;
  
  for(int i=0;i<tatom;i++){
    blist[i].clear();
  }
  delete[] blist;
  blist=NULL;



  if(check==NULL)check = new bool[gx.size()];
  for(int i=0;i<gx.size();i++){
    check[i]=true;
  }

  atomnum = gx.size();

}

void Grained::rnd_cut(){
  double tmp1 = (1.61+0.2*rnd());
  double tmp2 = (0.65+0.6*rnd());
  //double tmp2 = (0.9+0.3*rnd());
  ibondcut2=tmp1*tmp1; //bond length
  addcut2=tmp2*tmp2;
}

void Grained::BuildClist(){
  if(clist==NULL){
    clist = new vector<int>*[c_x];
    for(int i =0;i<c_x;i++){
      clist[i]=new vector<int>[c_y];
    }
  }
  else{
    for(int i=0;i<c_x;i++){
      for(int j=0;j<c_y;j++){
	clist[i][j].clear();
      }
    }
  }

  if(chead==NULL){
    chead = new double**[c_x];
    for(int i=0;i<c_x;i++){
      chead[i]=new double*[c_y];
      for(int j=0;j<c_y;j++){
	chead[i][j]=new double[4];
      }
    }
  }

  double lx = x_max-x_min;
  double ly = y_max-y_min;
  double dlx = lx/c_x;
  double dly = ly/c_y;

  for(int i=0;i<c_x;i++){
    for(int j=0;j<c_y;j++){
      chead[i][j][0]=x_min+dlx*i;
      chead[i][j][1]=x_min+dlx*i+dlx;
      chead[i][j][2]=y_min+dly*j;
      chead[i][j][3]=y_min+dly*j+dly;
    }
  }

  for(int i=0;i<tatom;i++){
    int cell_i=0;
    int cell_j=0;
    for(int ci=0;ci<c_x;ci++){
      for(int cj=0;cj<c_y;cj++){
	if(gx[i]>chead[ci][cj][0] && gx[i]<=chead[ci][cj][1]){
	  if(gy[i]>chead[ci][cj][2] && gy[i]<=chead[ci][cj][3]){
	    cell_i = ci;
	    cell_j = cj;
	  }
	}
      }
    }
    clist[cell_i][cell_j].push_back(i);
  }
  
  int acount=0;
  for(int ci=0;ci<c_x;ci++){
    for(int cj=0;cj<c_y;cj++){
      acount+=clist[ci][cj].size();
    }
  }
  

  if(myrank==0){
    cout << "Built Cell List!! \n";
    cout <<"Total atom number is "<<tatom<<" total atoms number of cells " <<acount<<endl;
  }

}

void Grained::SaveCoord(){
  sav_anum=atomnum;
  if(sav_x==NULL){
    sav_x=new double[sav_anum];
    sav_y=new double[sav_anum];

  }else{
    delete[] sav_x;
    delete[] sav_y;
    sav_x=new double[sav_anum];
    sav_y=new double[sav_anum];

  }
  
  int icount=0;
  for(int i=0;i<gx.size();i++){
    if(check[i]){
      
      sav_x[icount] = gx[i];
      sav_y[icount] = gy[i];
      icount++;
    }
  }

  if(icount!=sav_anum) cout <<"# of atoms does not match!\n";

}

void Grained::SetCoord(int nlocal, double **r){
  int lmpatom=nlocal;
  int *msize=new int[size];
  for(int i=0;i<size;i++) msize[i]=0;

  msize[myrank]=nlocal;
  
  local_isum(&lmpatom,1);
  local_isum(msize,size);
  
  if(myrank==0){printf("LAMMPS #of atoms %d\n",lmpatom); }

  sav_anum =lmpatom;
  tatom = lmpatom;
  atomnum = lmpatom;

  double *tmpx,*tmpy;

  tmpx=new double[sav_anum];
  tmpy=new double[sav_anum];

  if(sav_x==NULL){
    sav_x=new double[sav_anum];
    sav_y=new double[sav_anum];
  }else{
    delete[] sav_x;
    delete[] sav_y;
    sav_x=new double[sav_anum];
    sav_y=new double[sav_anum];
  }

  for(int i=0;i<sav_anum;i++){
    sav_x[i]=0.0;
    sav_y[i]=0.0;
    tmpx[i]=0.0;
    tmpy[i]=0.0;
  }

  int istart=0;
  for(int i=0;i<myrank;i++)istart+=msize[i];
  int iend = istart+nlocal;

  local_barrier();
  printf(" Myrank: %d, my start: %d, my end: %d\n",myrank,istart,iend);

  for(int i=istart;i<iend;i++){
    tmpx[i]=r[i-istart][0];
    tmpy[i]=r[i-istart][1];
  }
  
  local_dsum(tmpx,sav_anum);
  local_dsum(tmpy,sav_anum);

  for(int i=0;i<sav_anum;i++){
    sav_x[i]=tmpx[i];
    sav_y[i]=tmpy[i];

    //if(myrank==0)cout<<sav_x[i]<<" " <<sav_y[i]<<endl; 
  }

  delete[] msize;
  delete[] tmpx,tmpy;

  if(myrank==0)cout<<"Complete saving atoms from LAMMPS!\n";
  
}
