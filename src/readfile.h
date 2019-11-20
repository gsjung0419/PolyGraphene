//2013.10.07 G.S. Jung LAMM
#ifndef READ_FILE_H
#define READ_FILE_H
#include "myheader.h"
#include "atominfo.h"
class ReadFile
{
protected:
  AtomInfo *atominfo;
  void SetAtomInfo();

public:
  int blockx,blocky,blockz;
  ReadFile(){atominfo=NULL;
    blockx=1;blocky=1;blockz=1;
  }
  ~ReadFile(){if(atominfo!=NULL)delete atominfo;}

  //ReadCoord : The format is dump xyz (VMD Input format), timestep
  void ReadCoord(string _filename,int _tstep);

  //ReadPair : The format is dump custom x y z xx yy zz xy yz xz, number of cells
  void ReadCrack(string _filename,int _ncell);

  //xlo xhi..Format After readindg coordnates
  void ReadBox(string _filename);

  AtomInfo* GetAtomInfo(){return atominfo;}

  void RemoveTriangle(double width, double length);
  void RemoveYTriangle(double width, double length);
  void RemoveCenterx(double ycrackl, double width);
  void Translate(double x, double y, double z);

  void ReadPDB(string s);

};



#endif
