//2013.09.28 G.S. Jung @LAMM
#ifndef ATOM_INFO_H
#define ATOM_INFO_H
#include "myheader.h"

class AtomInfo
{
protected:

  
public:
  AtomInfo();
  ~AtomInfo();
  int numAtom;
  int numType;
  int blockx,blocky,blockz;

  double *box;
  double *subox;

  int *atype;
  int *ids;
  double *charge;

  double *x;
  double *y;
  double *z;

  double **r;
  double **vr;

  double *vx;
  double *vy;
  double *vz;

  void SetupArray();

  void SetAtomNum(int _numAtom);
  void SetAtomIds(int *_ids);
  void ReAtomNum(int _numAtom){numAtom=_numAtom;}
  void SetCoord(double *_x, double* _y,double* _z);
  void SetAtype(int *_atype);
  void SetBox(double * _h);

  int GetNumAtom(){return numAtom;}
  int GetNumType(){return numType;}
  int *GetAtype(){return atype;}
  int GetAtype(int i){return atype[i];}
  int *GetAtomIds(){return ids;}
  int GetAtomIds(int i){return ids[i];}
  double *GetCharge(){return charge;}
  double GetCharge(int i){return charge[i];}

  double * GetBox(){return box;}

  double *GetX(){return x;}
  double GetX(int i){return x[i];}
  double *GetY(){return y;}
  double GetY(int i){return y[i];}
  double *GetZ(){return z;}
  double GetZ(int i){return z[i];}
  double *GetVX(){return vx;}
  double GetVX(int i){return vx[i];}
  double *GetVY(){return vy;}
  double GetVY(int i){return vy[i];}
  double *GetVZ(){return vz;}
  double GetVZ(int i){return vz[i];}

  void AssignCell(std::vector<int> &clist, AtomInfo *gatominfo);
  void GenVMD(string s,string sa);

};

#endif
