//2013.10.15 G.S. Jung
#ifndef LAMMPS_COMMAND_H
#define LAMMPS_COMMAND_H
#include "myheader.h"
#include "atominfo.h"
#include "cellinfo.h"
#include "assemble.h"
#include "mpi.h"
#include "mympi.h"
#include "lammps.h"
#include "input.h"
#define CHNUM 200

using namespace LAMMPS_NS;
class LMPCommand
{
protected:
  LAMMPS *lammps;
  double timestep;
  int fix_num;
  int dump_num;
  vector <string> fix_command;
  vector <string> dump_command;

public:
  LMPCommand();
  ~LMPCommand();

  void Boundary(string s);
  void Compute(string s);
  void Dump(string s);
  void DumpModi(string s);
  void EAMInit(LAMMPS *_lammps);
  void REAXInit(LAMMPS *_lammps);
  void Fix(string s);
  void Minimize();
  void MinimizeISO();
  void MinimizeAniso(string cs);
  void ReadDataEAM(string s1,string s2);
  void ReadDataREAX(string s1,string s2);
  void ResetTime();
  void Run(int step);
  void Thermos(int tstep,string s);
  void TimeNeighbor(double f,double bin);
  void Units(string s);
  void UnFix(string s);
  void UnDump(string s);
  void Variable(string vname,string s);
  void Xdeform(int step,double erate);
};
#endif
