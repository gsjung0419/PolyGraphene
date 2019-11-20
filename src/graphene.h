//2013.10.03 G.S. Jung @LAMM
#ifndef GENERATE_GRAPHENE_H
#define GENERATE_GRAPHENE_H
#include "myheader.h"
#include "atominfo.h"
#include "gensolid.h"

class Graphene:public GenSolid{
protected:
  
public:
  Graphene();
  ~Graphene();
  AtomInfo* GenByCell(int _lat_x, int _lat_y, int _lat_z);
  AtomInfo* GenByDimension(double lx, double ly, double lz);
  void GenGB(double theta);
  void SetBondSheet(double _bondl, double _sheetz);
  void RemoveCenterx(double ycrackl, double width);
};

#endif
