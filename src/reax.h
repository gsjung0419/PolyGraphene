//2013.10.15 G.S. Jung
#ifndef LAMMPS_REAX_H
#define LAMMPS_REAX_H

#include "myheader.h"
#include <mpi.h>
#include "atominfo.h"
#include "readfile.h"
#include "lmpinput.h"
#include "gensolid.h"
#include "graphene.h"
#include "fcc.h"
#include "cellinfo.h"
#include "assemble.h"
#include "lammps.h"
#include "input.h"
#include "mympi.h"
#include "atom.h"
#include "domain.h"
#include "cell.h"
#include "library.h"
#include "lmpcommand.h"
#include "modify.h"
#include "compute.h"
#include "update.h"
#include "output.h"
#include "thermo.h"
#include "gslmp.h"

using namespace LAMMPS_NS;
class gsREAX:public gsLMP{
protected:

public:
  gsREAX();
  ~gsREAX();
  string* atoms;
  int ntype;

  void Init(int argc, char *argv[]);
  void GenGraphene(double _bodnl, double _sheetz, int _cx, int _cy, int _cz);
  void SetAtoms(string _filename);
  
  //From Virtual Functions
  void ReadDat();
  void ThermoUnit();
  void SetNPT(string cs);
  void ReadRestart(string _readfile);
};
#endif
