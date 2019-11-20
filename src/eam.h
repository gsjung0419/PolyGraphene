//2013.11.04 G.S. Jung
#ifndef LAMMPS_EAM_H
#define LAMMPS_EAM_H

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
class gsEAM:public gsLMP{
protected:

public:
  gsEAM();
  ~gsEAM();
  
  void Init(int argc, char *argv[]);
  void GenFCC(double _lc,double _rcut,int _cx, int _cy, int _cz);
  void GenFCC110(double _lc,double _rcut,int _cx, int _cy, int _cz);

  void SetMetal(string _aname,string _filename);
  // void CheckCoord();
  void Run(int tostep,int ostep);

  //From the Viratual Functions
  void ReadDat();
  void ThermoUnit();
  void SetNPT(string cs);
  void ReadRestart(string readfile);

};
#endif
