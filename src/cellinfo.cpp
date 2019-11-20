//2013.10.06 G.S. Jung @LAMM
#include "cellinfo.h"

CellInfo::CellInfo(){
  numAtom=0;
  rep_aid=-1;
  cell = NULL;
  atominfo =NULL;
  cen_x=0.0;
  cen_y=0.0;
  cen_z=0.0;
}

CellInfo::~CellInfo(){
  delete cell;
  delete atominfo;
  
}

void CellInfo::SetCell(double *_cbox){
  cell=new double[6];
  for(int i=0;i<6; i++){
    cell[i]=_cbox[i];
  }
  
  cen_x = 0.5*(cell[0]+cell[1]);
  cen_y = 0.5*(cell[2]+cell[3]);
  cen_z = 0.5*(cell[4]+cell[5]);


}

void CellInfo::AssignAtoms(AtomInfo *gatominfo){
  if(cell == NULL) {
    cout <<"Before You Assign Atoms, Set Cell Parameter! \n";
    exit(1);
  }
  
  int tNumAtom = gatominfo->GetNumAtom();
  double *x,*y,*z,*charge;
  int *atype,*ids;
  x = gatominfo->GetX();
  y = gatominfo->GetY();
  z = gatominfo->GetZ();

  vector <int> clist;
  for(int i=0;i<tNumAtom;i++){
    if(cell[0]<=x[i] && cell[1]>x[i]){
      if(cell[2]<=y[i] && cell[3]>y[i]){
	if(cell[4]<=z[i] && cell[5]>z[i]){
	  clist.push_back(i);
	}
      }
    }
  }  
  
  
  numAtom=clist.size();
  if(numAtom==0) {printf("There is a cell which doesn't have any atom num\n");
    exit(1);}
  

  atominfo=new AtomInfo;
  atominfo->SetAtomNum(numAtom);
  atominfo->AssignCell(clist, gatominfo);
  atominfo->SetBox(cell);
}

void CellInfo::SetId(int i_id){
  int *ids = atominfo->GetAtomIds();
  int numAtom = atominfo->GetNumAtom();
  for(int i=0;i<numAtom;i++){
    ids[i] = i_id+i;
  }
}
