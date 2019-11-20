//2013.10.03 G.S. Jung @LAMM
#ifndef GENERATE_FCC_H
#define GENERATE_FCC_H
#include "myheader.h"
#include "atominfo.h"
#include "gensolid.h"

class FCC:public GenSolid{
protected:

public:
  FCC();
  ~FCC();
  AtomInfo* GenByCell(int _lat_x, int _lat_y, int _lat_z);
  AtomInfo* GenByCell110(int _lat_x, int _lat_y, int _lat_z);
  AtomInfo* GenByCell111(int _lat_x, int _lat_y, int _lat_z);
  AtomInfo* GenByDimension(double lx, double ly, double lz);
  void SetLattice(double _lattice);
  void RemoveCenterx(double ycrackl, double width);
  void RemoveCenterxc(int x_width,int y_length);


};

#endif
