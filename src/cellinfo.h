//2013.10.06 G.S. Jung @LAMM
#ifndef CELL_INFORMATION_H
#define CELL_INFORMATION_H
#include "myheader.h"
#include "atominfo.h"

class CellInfo{
 protected:

 public:
  int numAtom;
  int rep_aid;
  double *cell;
  AtomInfo *atominfo;
  double cen_x,cen_y,cen_z;
  CellInfo();
  ~CellInfo();
  void SetCell(double *_cbox);
  void AssignAtoms(AtomInfo *gatominfo);
  int GetNumAtom(){return numAtom;}
  AtomInfo *GetAtomInfo(){return atominfo;}
  void SetId(int i_id);
  double GetX(){return cen_x;}
  double GetY(){return cen_y;}
  double GetZ(){return cen_z;}
  void SetRepAtom(int _id){rep_aid = _id;}
  int GetRepAtom(){return rep_aid;}
};

#endif 
