//2013.10.07 G.S. Jung @LAMM
#ifndef CELL_ASSEMBLE_H
#define CELL_ASSEMBLE_H
#include "myheader.h"
#include "atominfo.h"
#include "cellinfo.h"

class Assemble{
 protected:
  int numPart;
  int numCell;
  int numPair;
  int cx,cy,cz;
  int *cAtomnList;
  double *gbox;
  CellInfo **cellassemble;
  AtomInfo **atomassemble;
  int **pairs;
  int **cellpair;
  AtomInfo *gatominfo;

 public:
  Assemble();
  ~Assemble();
  void SetBox(double * _box);
  void SetCellNum(int _cx, int _cy, int _cz);
  void SetGlobalAtomInfo(AtomInfo *_gatominfo);
  void GenerateCell();
  CellInfo **GetCellInfo();
  void SetRepAtom();
  void SetPairCell();
  int GetNumCell(){return numCell;}
  int GetNumPair(){return numPair;}
  int**GetPairs(){return pairs;}
};

#endif 
