//2013.10.03 G.S. Jung @LAMM
#ifndef GRAIN_BOUNDARY_H
#define GRAIN_BOUnDARY_H
#include "myheader.h"
#include "mympi.h"
#include "atominfo.h"
#include "voro++.hh"
#include "graphene.h"
#include "readfile.h"


class Grained{
protected:
  
public:
  
  Grained();
  ~Grained();

  ofstream hfp;
  string hfname;
  int hstep;

  AtomInfo *atominfo;
  int myrank;
  int size;

  int particles;
  int n_x,n_y,n_z;
  int c_x,c_y;
  int atomnum;
  int tatom;
  int maxhanum;
  
  int sav_anum;
  double *sav_x,*sav_y;

  double *rotate;
  double ibondcut2;
  double addcut2;
  double mover_x;
  
  double x_min,x_max;
  double y_min,y_max;
  double z_min,z_max;

  vector<double> gx;
  vector<double> gy;
  vector<double> gz;

  vector<int> *blist;
  vector<int> **clist;
  double ***chead;

  bool *check;

  double rnd() {
    return double(rand())/RAND_MAX;
  }

  void GenVoro();  
  void GenVoroRe();
  void CheckBlist();
  void BuildBlist();
  void BuildBlist2();
  void BuildBlist3();
  void BuildClist();
  void GenLMP();
  void Readxyz1();
  void Readxyz2();
  void rnd_cut();

  void SaveCoord();
  void ReadCoord();
  void SetCoord(int nlocal, double **r);


  voro::container *con;

};

#endif
