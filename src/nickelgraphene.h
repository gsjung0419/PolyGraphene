//2013.10.03 G.S. Jung @LAMM
#ifndef NICKEL_GRAPHENE_H
#define NICKEL_GRAPHENE_H
#include "myheader.h"
#include "atominfo.h"
#include "gensolid.h"

class NickelGraphene:{
protected:


public:
  double *nickelbox;
  double *graphenebox;
  double *totalbox;

  NickelGraphene();
  ~NickelGraphene();

};

#endif
