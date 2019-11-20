//2013.10.06 G.S. Jung @LAMM
#ifndef GENERATE_SOLID_H
#define GENERATE_SOLID_H
#include "myheader.h"
#include "atominfo.h"

class GenSolid{
protected:
  AtomInfo *atominfo;
  double lattice;
  double mbondl;
  double bondl;
  double sheetz;
  int lat_x,lat_y,lat_z;
public:
  AtomInfo* GenFCC();
  AtomInfo* GenDIA();

  void SetBox(); //Only Orthogonal
  void RemoveBox(double *rbox);
  void RemovePoint(double *r);
  void CellDivision(int x, int y, int z);
  void MakeSpace(double dx, double dy, double dz);
  void RemoveTriangle(double width, double length);
  void RemoveYTriangle(double width, double length);
  
  //virtual AtomInfo* GenByCell(int _lat_x, int _lat_y, int _lat_z)=0;
  //virtual AtomInfo* GenByDimension(double lx, double l_y, double l_z)=0;
  //virtual void RemoveCenterx(double ycrackl, double width);

};

#endif
