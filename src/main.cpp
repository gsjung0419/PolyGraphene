//2013.11.04 G.S. Jung
#include "myheader.h"
#include <mpi.h>
#include "atominfo.h"
#include "readfile.h"
#include "gensolid.h"
#include "graphene.h"
#include "fcc.h"
#include "assemble.h"
#include "input.h"
#include "mympi.h"
#include "cell.h"
#include "modify.h"
#include "eam.h"
#include "reax.h"
#include "analysis.h"
#include "voro++.hh"
#include "grained.h"


double rnd() {return double(rand())/RAND_MAX;}
const int particles=5;

MPI_Comm myComm;
using namespace LAMMPS_NS;

void GenVoro();

int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  myComm=MPI_COMM_WORLD;
  
  ReadFile *readf;
  AtomInfo *atominfo;
  LMPInput *lmpinp;
  Grained *gb;
  
  gb = new Grained;
  gb->GenVoro();

  using namespace LAMMPS_NS;  
  LAMMPS *lmp = new LAMMPS(argc,argv,MPI_COMM_WORLD);
  for(int i=0;i<10000;i++){
   if(i!=0){
     gb->rnd_cut();
     gb->Readxyz1();//gen from voro
     gb->GenLMP();
   }

   delete lmp;
   lmp = new LAMMPS(argc,argv,MPI_COMM_WORLD);
   lmp->input->file("relax.in");
   
    gb->rnd_cut();
    gb->Readxyz2();//add voro
  }
  
  gb->Readxyz1();//gen from voro

  MPI_Finalize();
  return 0;

}

