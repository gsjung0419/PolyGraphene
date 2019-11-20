//2013.11.04 G.S. Jung
#ifndef LAMMPS_GSLMP_H
#define LAMMPS_GSLMP_H

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

using namespace LAMMPS_NS;
class gsLMP
{
protected:

  LMPInput *lmpinput;
  Thermo *thermo;
  FCC *fcc;
  Graphene *graphene;

  AtomInfo *atominfo;
  
  bool on_crack;
  bool on_erate;

  ofstream ofstress;
  ofstream ofthermo;
  ofstream oftslaw;
  int myrank;

  double etotal0;
  double lx0,ly0,lz0;

  void LMPVMD();
  void LMPCFG();
  void TSLAW();

  void SetComputes();
  void Deform(string cs);

 public:
  LMPCommand *cmd;
  LAMMPS *lammps;
  Cells *cells;



  string potential;
  string metal;
  string restartfile;
  string thermoout;
  string stressout;
  string vmdout;
  string cfgout;
  string tslaw;

  int deformStep;
  int thermoStep;
  int rescaleStep;
  int outputStep;
  int vmdStep;
  int cfgStep;
  int scanStep;
  int tslawStep;
  int numcell;
  int mstep;

  bool on_vmd;
  bool on_cfg;
  bool on_deform;
  bool on_strain;
  bool on_stmin;
  bool on_tslaw;

  int cx,cy,cz;

  double lattice;
  double bondl;
  double sheetz;
  double timestep;
  double temperature;
  double pressure;
  double erate;
  double delx;
  double crackL;
  double rcut;
  double *boxhi;
  double *boxlo;

  //gsLMP();
  //~gsLMP();

  void Crack(int _width,int _length);

  void SetThermo(int _thermostep, double _temperature, double _pressure,string s);
  void Equilibrium(int tostep,double _temperature,double _pressure);


  void WriteRestart(string writefile);

  void SetVMD(int _ostep, string _s);
  void SetCFG(int _ostep, string _s);

  void SetNPTISO();

  void SetStMin(double perstrain, string _s);
  void SetStrain(double _erate,int _dstep, int _scan,string _s);
  //void SetDeform(double _erate,int _dstep, int _scan,string _s);
  void SetTSLaw(double _rcut,int numcell,int _tsnum, string s);

  void RunStMin(double estrain,string cs);
  //void RunDeform(double estrain,string cs);
  void RunStrain(double estrain,string cs);

  void PrintThermo(int _step);
  void PrintStress(int _step);
  void PrintAtomicStress(int _step);

  virtual void ReadDat(){};
  virtual void ThermoUnit(){};
  virtual void SetNPT(string cs){};
  virtual void ReadRestart(string readfile){};
};
#endif
