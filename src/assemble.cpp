//2013.10.10 G.S. Jung @LAMM
#include "assemble.h"

Assemble::Assemble()
{
  numCell=0;
  numPart=0;
  numPair=0;
  cx=0;cy=0;cz=0;
  cellassemble=NULL;
  atomassemble=NULL;
  gbox=new double[6];
  pairs=NULL;
  gatominfo=NULL;
  cAtomnList=NULL;
}

Assemble::~Assemble()
{
  for(int i=0;i<numCell;i++) delete cellassemble[i];
  for(int i=0;i<numPart;i++) delete atomassemble[i];
  for(int i=0;i<numPair;i++) delete pairs[i];


  delete[] cAtomnList;
  delete[] gbox;

  delete [] cellassemble;
  delete [] atomassemble;
  delete [] pairs;
  //delete gatominfo; //You Cannot Delete. 
}

void Assemble::SetCellNum(int _cx, int _cy, int _cz){
  cx = _cx;
  cy = _cy;
  cz = _cz;
  numCell = cx*cy*cz;


}

void Assemble::SetGlobalAtomInfo(AtomInfo *_gatominfo){
  gatominfo = _gatominfo;
  double *h = gatominfo->GetBox();
  for(int i=0;i<6;i++){gbox[i] = h[i];}

}

void Assemble::GenerateCell(){
  if(cx==0 || cy==0 || cz==0){
    cout <<"Error: SetCellNum before GenerateCell\n";
    exit(1);
  }
  if(gatominfo==NULL){
    cout <<"Error: SetGlobalAtomInfo before GenerateCell\n";
    exit(1);
  }

  //cout <<"before cellassemble" <<endl;  
  cellassemble = new CellInfo*[numCell];
  for(int i=0;i<numCell;i++) cellassemble[i]=new CellInfo;

  double lx = (gbox[1]-gbox[0])/cx;
  double ly = (gbox[3]-gbox[2])/cy;
  double lz = (gbox[5]-gbox[4])/cz;
  
  cAtomnList = new int[numCell]; //AtomNumberList

  //cout <<"before Assignatoms" <<endl;  
  for(int i=0;i<cx;i++){
    for(int j=0;j<cy;j++){
      for(int k=0;k<cz;k++){
	int index = (i*cy*cz + j*cz + k);
	double cbox[6]={0.0};
	cbox[0]=gbox[0]+lx*i;
	cbox[1]=gbox[0]+lx*(i+1);
	cbox[2]=gbox[2]+ly*j;
	cbox[3]=gbox[2]+ly*(j+1);
	cbox[4]=gbox[4]+lz*k;
	cbox[5]=gbox[4]+lz*(k+1);

	cellassemble[index]->SetCell(cbox);
	cellassemble[index]->AssignAtoms(gatominfo);
	cAtomnList[index] = cellassemble[index]->GetNumAtom();	
      }
    }
  }
  
  int *iid = new int[numCell];
  iid[0] = 1;
  
  for(int i=1;i<numCell;i++){
    iid[i]=cAtomnList[i-1] + iid[i-1];
  }

  for(int i=0;i<numCell;i++){
    cellassemble[i]->SetId(iid[i]);
  }


  //cout <<"before setpaircell" <<endl;  
  if(cx==2){
    SetPairCell();
    SetRepAtom();
  }

  delete[] iid;
}

CellInfo **Assemble::GetCellInfo(){
  if(cellassemble==NULL){
    cout <<"GenerateCell Before GetCellAssemble\n";
    exit(1);
  }
  return cellassemble;

}

void Assemble::SetRepAtom(){
  double lx=0.5*(gbox[1]-gbox[0])+gbox[0];
  double ly=0.5*(gbox[3]-gbox[2])+gbox[2];
  double lz=0.5*(gbox[5]-gbox[4])+gbox[4];

  for(int i=0; i<numPair;i++){
    int icell = pairs[i][0];
    int jcell = pairs[i][1];
    
    AtomInfo *icellatoms = cellassemble[icell]->GetAtomInfo();
    AtomInfo *jcellatoms = cellassemble[jcell]->GetAtomInfo();
    
    int inumAtom = icellatoms->GetNumAtom();
    int jnumAtom = jcellatoms->GetNumAtom();
    
    int rep_ia=-1;
    int rep_ja=-1;
    double dist=-1;
    
    for(int ia=0;ia<inumAtom;ia++){
      for(int ja=0;ja<jnumAtom;ja++){
	double xi = icellatoms->GetX(ia);
	double xj = jcellatoms->GetX(ja);

	double yi = icellatoms->GetY(ia);
	double yj = jcellatoms->GetY(ja);

	double zi = icellatoms->GetZ(ia);
	double zj = jcellatoms->GetZ(ja);
	
	double l2 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) +(zi-zj)*(zi-zj);
	
	if(ia==0 && ja==0){
	  dist = l2;
	  rep_ia=ia;
	  rep_ia=ja;
	}
	else if(dist>l2){
	  dist=l2;
	  rep_ia=ia;
	  rep_ja=ja;
	}

      }
    }

    cellassemble[icell]->SetRepAtom(rep_ia);
    cellassemble[jcell]->SetRepAtom(rep_ja);

  }

  for(int i=0; i<numCell;i++){
    int rep_id = cellassemble[i]->GetRepAtom();
    AtomInfo* atominfo=cellassemble[i]->GetAtomInfo();
    int id = atominfo->GetAtomIds(rep_id);
    //cout << i << " " <<rep_id<< " " << atominfo->GetX(rep_id) << " " <<atominfo->GetY(rep_id)<<" "  <<atominfo->GetZ(rep_id) <<endl;
  }
}

void Assemble::SetPairCell(){
  double ly = (gbox[3] - gbox[2])/cy;
  pairs = new int*[cy];
  for(int i=0;i<cy;i++) pairs[i]=new int[2];

  int cpair=0;
  for(int i=0; i<numCell;i++){
    double pair1 = cellassemble[i]->GetY();
    
    for(int j=i+1; j<numCell;j++){
      double pair2 = cellassemble[j]->GetY();
      double diff = fabs(pair1-pair2);
      
      if(diff < 0.5*ly ){
	pairs[cpair][0]=i;
	pairs[cpair][1]=j;
	cpair++;
      }
    }
  }

  numPair=cpair;
  
  for(int i=0;i<cy;i++){
    //cout <<pairs[i][0] << " " <<pairs[i][1]<<endl;
  }

}

