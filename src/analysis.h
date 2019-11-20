//2013.10.03 G.S. Jung @LAMM
#ifndef MY_ANALYSIS_H
#define MY_ANALYSIS_H
#include "myheader.h"
#include "atominfo.h"
#include "gensolid.h"

class Analysis{
protected:

public:
  int numfile;
  int numline;
  int numdat;
  int skipline;
  string filehead;

  double **mydata;

  Analysis();
  ~Analysis();

  void NewData();
  void CloseData();
  void Average();
  void ThermoLog();

  

};

#endif
