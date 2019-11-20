//2013.09.30 G.S. Jung
#ifndef LAMMPS_INPUT_H
#define LAMMPS_INPUT_H
#include "myheader.h"
#include "atominfo.h"

#define EAM 0
#define REAX 1
#define CHARM 2
#define CVFF 3
#define EAM_REAX 4


class LMPInput
{
protected:
  AtomInfo *atominfo;

 public:					
  int snum;
  int numAtomType;
  int poType;

  LMPInput();
  ~LMPInput();
  void GenInput(string s);
  void GenDat(string s);
  void GenDatTypeNum(string s, int num);
  void SetAtomInfo(AtomInfo* _atominfo){
    atominfo=_atominfo;
  }
  void GenVMD(string s);
};
#endif
