//2013.10.19 G.S. Jung @LAMM
#ifndef CELL_H
#define CELL_H
#include "myheader.h"
#include "atominfo.h"
#include "atom.h"
#include "lammps.h"
#include "library.h"
#include "domain.h"
#include "input.h"
#include "compute.h"
#include "modify.h"


using namespace LAMMPS_NS;
class Cells{
 protected:
  ofstream fp;
  ofstream tsf;
  int cx,cy,cz;
  int natoms;
  int nlocal;
  int fp_count;

  double *boxlo;
  double *boxhi; 
  double *sublo;
  double *subhi;
  LAMMPS *lammps;
  Atom *atom; 

  bool seperation;
  double **dist0;

  //void UpdateGlobal();
  void UpdateCell();
  void UpdatePair();
  void UpdateCellAtom();
  
  void UpdateLocal();
  void UpdateLocalPair();
  void UpdateLocalAtom();

  int myrank;
  int framenum;
 public:
  int numCell;
  int numPair;

  int *gtags;
  int *ltags;
  int *rep_atom;
  int *prep_atom;
  int *listnumAtom;
  int *lpair;
  int *nnum;

  int **pairs;
  int **atoms;
  int **nlist;

  double rcut;
  double crackL;
  double **gatoms;
  double **cbox;
  double **pdist;
  double **sdist;
  double **repstress;
  double *avolume;
  int *dcount;

  void InitCell(int _cx,int _cy,int _cz,LAMMPS *_lammps);
  //void UpdateRepAtom();
  void UpdateLocalRepAtom();
  void UpdateLocalRepStress();
  void UpdateGroup(string s);
  void UnGroup(string s);
  void AnalyseTSLaw(string s);
  void AvgTSLaw(int total, int numavg, string s);
  void SetRcutCrackL(double _rcut){
    rcut = _rcut;
  }
  
  Cells();
  ~Cells();

};

#endif 
