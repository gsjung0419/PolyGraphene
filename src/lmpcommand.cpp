//2013.11.04 G.S. Jung @LAMM
#include "lmpcommand.h"

LMPCommand::LMPCommand(){
  lammps=NULL;
  fix_num=5000;			       
  dump_num=2000;
  timestep=0.0;
}

LMPCommand::~LMPCommand(){
  lammps=NULL;
}

void LMPCommand::Boundary(string s){
  char com[CHNUM]={NULL};
  sprintf(com,"change_box  \t all boundary %s",s.c_str());
  lammps->input->one(com);
}

void LMPCommand::Compute(string s){
  char command1[CHNUM]={NULL};
  sprintf(command1,"compute     \t %s",s.c_str());

  lammps->input->one(command1);
}

void LMPCommand::Dump(string s){
  dump_command.push_back(s);
  dump_num++;
  char command1[CHNUM]={NULL};
  sprintf(command1,"dump        \t %d %s",dump_num,s.c_str());
  lammps->input->one(command1);
}

void LMPCommand::DumpModi(string s){
  char command1[CHNUM]={NULL};
  sprintf(command1,"dump_modify \t %d %s",dump_num,s.c_str());
  lammps->input->one(command1);
}

void LMPCommand::EAMInit(LAMMPS *_lammps){
  //Initialization for EAM Potential
  lammps=_lammps;
  Units("metal");
  lammps->input->one("dimension   \t 3");
  lammps->input->one("boundary    \t p p p");
  lammps->input->one("atom_style  \t atomic");
}

void LMPCommand::REAXInit(LAMMPS *_lammps){
  //Initialization for Reactive ForceField
  lammps=_lammps;
  Units("real");
  lammps->input->one("dimension   \t 3");
  lammps->input->one("boundary    \t p p p");
  lammps->input->one("atom_style  \t charge");
}

void LMPCommand::Fix(string s){
  char command1[CHNUM]={NULL};
  sprintf(command1,"fix         \t %d %s",fix_num,s.c_str());
  lammps->input->one(command1);
  fix_command.push_back(s);
  fix_num++;
}

void LMPCommand::Minimize(){

  lammps->input->one("min_style   \t cg");
  lammps->input->one("minimize    \t 1.0e-10 1.0e-10 100000 100000");
}

void LMPCommand::MinimizeISO(){
  lammps->input->one("fix         \t 100 all box/relax iso 0.0");
  lammps->input->one("min_style   \t cg");
  lammps->input->one("minimize    \t 1.0e-10 1.0e-10 100000 100000");
  lammps->input->one("unfix       \t 100");
}

void LMPCommand::MinimizeAniso(string cs){
  string tmpx="x";
  string tmpy="y";
  string tmpz="z";
  
  if(cs==tmpx){
    lammps->input->one("fix         \t 100 all box/relax y 0 z 0");
    lammps->input->one("min_style   \t cg");
    lammps->input->one("minimize    \t 1.0e-10 1.0e-10 100000 100000");
    lammps->input->one("unfix       \t 100");
  }

  if(cs==tmpy){
    lammps->input->one("fix         \t 100 all box/relax x 0  z 0");
    lammps->input->one("min_style   \t cg");
    lammps->input->one("minimize    \t 1.0e-10 1.0e-10 100000 100000");
    lammps->input->one("unfix       \t 100");
  }

  if(cs==tmpz){
    lammps->input->one("fix         \t 100 all box/relax x 0 y 0");
    lammps->input->one("min_style   \t cg");
    lammps->input->one("minimize    \t 1.0e-10 1.0e-10 100000 100000");
    lammps->input->one("unfix       \t 100");
  }

}

void LMPCommand::Run(int step){
  char command1[CHNUM]={NULL};
  sprintf(command1,"run         \t %d",step);
  lammps->input->one(command1);
}

void LMPCommand::ReadDataEAM(string s1,string s2){
  char command1[CHNUM]={NULL};
  char command2[CHNUM]={NULL};

  sprintf(command1,"read_data   \t %s",s1.c_str());
  sprintf(command2,"pair_coeff  \t * * %s",s2.c_str());
  lammps->input->one(command1);
  lammps->input->one("pair_style  \t eam/alloy");
  lammps->input->one(command2);

}

void LMPCommand::ReadDataREAX(string s1,string s2){
  char command1[CHNUM]={NULL};
  char command2[CHNUM]={NULL};

  sprintf(command1,"read_data   \t %s",s1.c_str());
  sprintf(command2,"pair_coeff  \t * * %s",s2.c_str());
  lammps->input->one(command1);
  lammps->input->one("pair_style  \t reax");
  lammps->input->one(command2);

}

void LMPCommand::ResetTime(){
  lammps->input->one("reset_timestep 0");
}

void LMPCommand::Thermos(int tstep,string s){
  char command1[CHNUM]={NULL};
  char command2[CHNUM]={NULL};
  sprintf(command1,"thermo      \t %d",tstep);
  sprintf(command2,"thermo_style\t %s",s.c_str());
  lammps->input->one(command1);
  lammps->input->one(command2);
}

void LMPCommand::TimeNeighbor(double f, double bin){
  char command1[CHNUM]={NULL};
  char command2[CHNUM]={NULL};
  sprintf(command1,"timestep    \t %f",f);
  sprintf(command2,"neighbor    \t %f bin",bin);
  lammps->input->one(command1);
}

void LMPCommand::Units(string s)
{  char command[CHNUM]={NULL};
  sprintf(command,"units       \t %s",s.c_str());
  lammps->input->one(command);
}

void LMPCommand::UnFix(string s){
  int t_fix = fix_command.size();
  char command[CHNUM]={NULL};
  bool found=false;
  for(int i=0;i<t_fix;i++){
    if(fix_command[i] == s){
      int find_num=fix_num-t_fix+i;
      sprintf(command,"unfix       \t %d",find_num);      
      lammps->input->one(command);
      fix_command[i]="N";
      found =true;
    }
  }

  if(!found){
    cout <<"There is no such fix command: "<<s<<endl;
  }

}

void LMPCommand::UnDump(string s){
  int t_dump = dump_command.size();
  char command[CHNUM]={NULL};
  bool found=false;
  for(int i=0;i<t_dump;i++){
    if(dump_command[i] == s){
      int find_num=dump_num-t_dump+i+1;
      sprintf(command,"undump      \t %d",find_num);      
      lammps->input->one(command);
      dump_command[i]="N";
      found =true;
    }
  }

  if(!found){
    cout <<"There is no such dump command: "<<s<<endl;
  }

}

void LMPCommand::Variable(string vname,string s){
  char command[CHNUM]={NULL};
  sprintf(command,"variable    \t %s %s",vname.c_str(),s.c_str());
  lammps->input->one(command);
}

void LMPCommand::Xdeform(int step, double erate){
  char command[CHNUM]={NULL};
  sprintf(command,"fix         \t %d all deform %d x erate %f units box remap x",fix_num,step,erate);
  lammps->input->one(command);
}


