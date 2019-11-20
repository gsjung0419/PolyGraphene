//2013.10.19 G.S. Jung @LAMM
#include "cell.h"
#include "mympi.h"
Cells::Cells(){
  numCell=0;
  numPair=0;
  cx=1;
  cy=1;
  cz=1;
  framenum=0;
  fp_count=0;

  crackL=0.0;
  rcut=0.0;

  lpair=NULL;

  boxlo=NULL;
  boxhi=NULL;
  sublo=NULL;
  subhi=NULL;
  lammps=NULL;
  atom=NULL;
  dcount=NULL;

  rep_atom=NULL;
  prep_atom=NULL;
  pairs=NULL;
  atoms=NULL;
  listnumAtom=NULL;
  gtags=NULL;
  
  cbox=NULL;
  pdist=NULL;
  sdist=NULL;
  gatoms=NULL;
  
  nlist=NULL;
  nnum=NULL;

  repstress=NULL;
  avolume=NULL;
  seperation=true;


}

Cells::~Cells(){
  delete[] rep_atom;
  delete[] listnumAtom;
  for(int i=0;i<numPair;i++) delete []pairs[i];
  if(atoms!=NULL) for(int i=0;i<numCell;i++) if(atoms[i]!=NULL)delete []atoms[i];
  if(cbox!=NULL)for(int i=0;i<numCell;i++) delete []cbox[i];
  if(pdist!=NULL)for(int i=0;i<numCell;i++) delete []pdist[i];
  if(pairs!=NULL)delete [] pairs;
  if(atoms!=NULL)delete [] atoms;
  if(cbox!=NULL)delete [] cbox;
  if(pdist!=NULL)delete [] pdist;
  if(myrank==0)fp.close();
  if(myrank==0)tsf.close();
}

void Cells::InitCell(int _cx,int _cy,int _cz,LAMMPS *_lammps){
  cx=_cx;
  cy=_cy;
  cz=_cz;
  
  lammps = _lammps;

  atom=lammps->atom;
  boxlo=lammps->domain->boxlo;
  boxhi=lammps->domain->boxhi;
  sublo=lammps->domain->boxlo;
  subhi=lammps->domain->boxhi;

  for(int i=0;i<atom->natoms;i++){
    //cout << atom->tag[i]<< " "<<atom->x[i][0] << " " << atom->x[i][1]<< " " << atom->x[i][2]<<endl;
  }

  natoms = lammps_get_natoms(lammps);
  nlocal = lammps->atom->nlocal;

  dcount = new int[1000];
  for(int i=0;i<1000;i++)dcount[i]=0;

  myrank = get_global_rank();

  if(myrank==0) fp.open("id_history.dat");
  if(myrank==0) tsf.open("tslaw.dat");
  if(myrank==0) tsf <<"cell#1 " << "yposition2 " << "distx(A)3 " << "separation(A)4 "<< "Avg_atomicVolume(A^3)5 " <<"hydrostress(GPa)6 "<< "pxx(GPa)7 " <<"pyy(GPa)8 "<<"pzz(GPa)9 "<< "vonmisses(GPa)10 "<<endl;

}

void Cells::UpdateLocal(){
  if(gatoms==NULL){
    gatoms=new double*[nlocal]; for(int i=0;i<nlocal;i++) gatoms[i]=new double[3];
    gtags=new int[nlocal];

    boxlo=new double[3];
    boxhi=new double[3];
    sublo=new double[3];
    subhi=new double[3];
    
  }else{
    for(int i=0;i<nlocal;i++){delete[] gatoms[i];}
    delete[] gatoms;
    delete[] gtags;
    delete[] boxlo;
    delete[] boxhi;
    delete[] sublo;
    delete[] subhi;

    nlocal = lammps->atom->nlocal;
    gatoms=new double*[nlocal]; for(int i=0;i<nlocal;i++) gatoms[i]=new double[3];
    gtags=new int[nlocal];

    boxlo=new double[3];
    boxhi=new double[3];
    sublo=new double[3];
    subhi=new double[3];
  }

  atom=lammps->atom;

  for(int i=0;i<3;i++){
    boxlo[i]=lammps->domain->boxlo[i];
    boxhi[i]=lammps->domain->boxhi[i];
    sublo[i]=lammps->domain->sublo[i];
    subhi[i]=lammps->domain->subhi[i];
  }

  for(int i=0;i<nlocal;i++){
    //cout << i<< " "<<atomg[3*i] << " " << atomg[3*i+1]<< " " << atomg[3*i+2]<<endl;
    gatoms[i][0]=atom->x[i][0];
    gatoms[i][1]=atom->x[i][1];
    gatoms[i][2]=atom->x[i][2];
    gtags[i]=atom->tag[i];
    //cout <<i << " " <<gtags[i]<<" " <<atom->tag[i]<<endl;
  }

}

void Cells::UpdateLocalPair(){
  if(lpair==NULL){
    lpair = new int[numPair];
    for(int i=0;i<numPair;i++) lpair[i]=0;
  }else{
    delete[] lpair;
    lpair = new int[numPair];
    for(int i=0;i<numPair;i++) lpair[i]=0;
  }

  int *ltmp =new int[numPair];
  for(int i=0;i<numPair;i++) ltmp[i]=0;

  for(int i=0;i<numPair;i++){
    double center[3]={0.0};
    int icell=pairs[i][0];
    
    center[0] = 0.5*(cbox[icell][0]+cbox[icell][1]);
    center[1] = 0.5*(cbox[icell][2]+cbox[icell][3]);
    center[2] = 0.5*(cbox[icell][4]+cbox[icell][5]);

    if(sublo[0]<center[0] && center[0]<subhi[0]){
      if(sublo[1]<center[1] && center[1]<subhi[1]){
	if(sublo[2]<center[2] && center[2]<subhi[2]){
	  lpair[i] = 1;
	}
      }
    }
  }

  global_isum(lpair,ltmp,numPair);
  
  bool check=true;
  for(int i=0;i<numPair;i++){
    if(ltmp[i]!=1){
      cout <<"Error : the local cell is not properly generated\n";
      for(int i=0;i<numPair;i++){
	cout << i << " "<< ltmp[i]<<endl;
      }
      exit(1);
    }
  }

  delete[] ltmp;
}

void Cells::UpdateLocalAtom(){
  if(atoms==NULL){
    atoms = new int*[numCell];
    for(int i=0;i<numCell;i++) atoms[i]=NULL;
    listnumAtom = new int[numCell];
  }else{
    for(int i=0;i<numCell;i++) if(atoms[i]!=NULL) delete[] atoms[i];
    delete [] listnumAtom;
    delete [] atoms;
    
    atoms=new int*[numCell];
    for(int i=0;i<numCell;i++) atoms[i]=NULL;
    listnumAtom = new int[numCell];
  }

  for(int i=0;i<numCell;i++){
    
    listnumAtom[i]=0;

    bool check=false;
    for(int j=0;j<numPair;j++){
      if(lpair[j] == 1) {
	int ic = pairs[j][0];
	int jc = pairs[j][1];
	if(ic == i || jc == i)check=true;
      }  
    }

    if(check){
      int count=0;
      for(int j=0;j<nlocal;j++){
	//if(cbox[i][0] <= atom->x[j][0] && cbox[i][1]> atom->x[j][0]){
	// if(cbox[i][2] <= atom->x[j][1] && cbox[i][3]> atom->x[j][1]){
	//  if(cbox[i][4] <= atom->x[j][2] && cbox[i][5]> atom->x[j][2]){
	if(cbox[i][0] <= gatoms[j][0] && cbox[i][1]> gatoms[j][0]){
	  if(cbox[i][2] <= gatoms[j][1] && cbox[i][3]> gatoms[j][1]){
	    if(cbox[i][4] <= gatoms[j][2] && cbox[i][5]> gatoms[j][2]){
	      count++;
	    }
	  }
	}
      }
      
      atoms[i]=new int[count];
      listnumAtom[i]=count;
      
      count=0;
      
      for(int j=0;j<nlocal;j++){
	//if(cbox[i][0] <= atom->x[j][0] && cbox[i][1]>atom->x[j][0]){
	//if(cbox[i][2] <= atom->x[j][1] && cbox[i][3]>atom->x[j][1]){
	//if(cbox[i][4] <=atom->x[j][2] && cbox[i][5]>atom->x[j][2]){
	if(cbox[i][0] <= gatoms[j][0] && cbox[i][1]> gatoms[j][0]){
	  if(cbox[i][2] <= gatoms[j][1] && cbox[i][3]> gatoms[j][1]){
	    if(cbox[i][4] <= gatoms[j][2] && cbox[i][5]> gatoms[j][2]){
	      atoms[i][count]=j;
	      count++;
	    }
	  }
	}
      }
    }
  }  
  
  int nrank = get_global_size();
  
  for (int i=0;i<numCell;i++){

    if(listnumAtom[i]==0 && nrank==1 ){ printf("Too many cell there are one or more cell which have no atoms!\n");
    cout<<i<<" " <<listnumAtom[i]<<endl;      
      exit(1);
    }
  }
}

void Cells::UpdateCell(){
  numCell=cx*cy*cz;
  //cout <<"Num Cell : "<<numCell;
  
  if(cbox==NULL){
    cbox=new double*[numCell];
    for(int i=0;i<numCell;i++){
      cbox[i]=new double[6];
    }
    
  }else{
    for(int i=0;i<numCell;i++) delete[] cbox[i];
    delete[] cbox;
    cbox=new double*[numCell];
      for(int i=0;i<numCell;i++){
	cbox[i]=new double[6];
      }
  }

  double lx=(boxhi[0]-boxlo[0])/cx;
  double ly=(boxhi[1]-boxlo[1])/cy;
  double lz=(boxhi[2]-boxlo[2])/cz;

  double lxo=boxlo[0];
  double lyo=boxlo[1];
  double lzo=boxlo[2];

  for(int i=0;i<cx;i++){
    for(int j=0;j<cy;j++){
      for(int k=0;k<cz;k++){
	int index = (i*cy*cz + j*cz + k);
	cbox[index][0]=lxo+lx*i;
	cbox[index][1]=lxo+lx*(i+1);
	cbox[index][2]=lyo+ly*j;
	cbox[index][3]=lyo+ly*(j+1);
	cbox[index][4]=lzo+lz*k;
	cbox[index][5]=lzo+lz*(k+1);
      }
    }
  }  

  //cout <<"Complete cbox gen\n";
}

void Cells::UpdatePair(){

  if(pairs==NULL){
    pairs = new int*[cy];
    for(int i=0;i<cy;i++) pairs[i]=new int[2];
  }else{
    for(int i=0;i<numPair;i++) delete pairs[i];
    delete pairs;
    pairs = new int*[cy];
    for(int i=0;i<cy;i++) pairs[i]=new int[2];
  }

  double ly = (boxhi[1]-boxlo[1])/cy;
  int cpair=0;
  for(int i=0; i<numCell;i++){
    double pair1 = 0.5*(cbox[i][2]+cbox[i][3]);
    for(int j=i+1; j<numCell;j++){
      double pair2 = 0.5*(cbox[j][2]+cbox[j][3]);
      double diff = fabs(pair1-pair2);
      if(diff < 0.5*ly){
	pairs[cpair][0]=i;
	pairs[cpair][1]=j;
	cpair++;
      }
    }
  }

  for(int i=0; i<numCell;i++){
    //printf("Cell %d : %f %f %f %f %f %f\n",i,cbox[i][0],cbox[i][1],cbox[i][2],cbox[i][3],cbox[i][4],cbox[i][5]);
  }


  numPair=cpair;
}

void Cells::UpdateCellAtom(){
  if(atoms==NULL){
    atoms = new int*[numCell];
    listnumAtom = new int[numCell];
  }else{
    for(int i=0;i<numCell;i++) delete[] atoms[i];
    delete [] listnumAtom;
    delete [] atoms;
    
    atoms=new int*[numCell];
    listnumAtom = new int[numCell];
  }

  for(int i=0;i<numCell;i++){
    int count=0;
    for(int j=0;j<natoms;j++){
      //if(cbox[i][0] <= atom->x[j][0] && cbox[i][1]> atom->x[j][0]){
      //if(cbox[i][2] <= atom->x[j][1] && cbox[i][3]> atom->x[j][1]){
      //if(cbox[i][4] <= atom->x[j][2] && cbox[i][5]> atom->x[j][2]){
      if(cbox[i][0] <= gatoms[j][0] && cbox[i][1]> gatoms[j][0]){
      if(cbox[i][2] <= gatoms[j][1] && cbox[i][3]> gatoms[j][1]){
        if(cbox[i][4] <= gatoms[j][2] && cbox[i][5]> gatoms[j][2]){
	    count++;
	  }
	}
      }
    }
      
    atoms[i]=new int[count];
    listnumAtom[i]=count;
    
    count=0;
    
    for(int j=0;j<natoms;j++){
      //if(cbox[i][0] <= atom->x[j][0] && cbox[i][1]>atom->x[j][0]){
      //if(cbox[i][2] <= atom->x[j][1] && cbox[i][3]>atom->x[j][1]){
      //if(cbox[i][4] <=atom->x[j][2] && cbox[i][5]>atom->x[j][2]){
	if(cbox[i][0] <= gatoms[j][0] && cbox[i][1]> gatoms[j][0]){
	if(cbox[i][2] <= gatoms[j][1] && cbox[i][3]> gatoms[j][1]){
	if(cbox[i][4] <= gatoms[j][2] && cbox[i][5]> gatoms[j][2]){
	    atoms[i][count]=j;
	    count++;
	  }
	}
      }
    }
  }  
  
  for (int i=0;i<numCell;i++){
    //cout<<i<<" " <<listnumAtom[i]<<endl;
    if(listnumAtom[i]==0){ printf("Too many cell there are one or more cell which have no atoms!");
      exit(1);
    }
  }
  
  

}

/*void Cells::UpdateRepAtom(){
  UpdateGlobal();
  //cout<<"Up Global\n";
  //for(int i=0;i<3;i++) cout<<boxhi[i]<<" "<<endl;    
  UpdateCell();
  //cout<<"Up Cell\n";
  UpdatePair();
  UpdateCellAtom();

  if(rep_atom==NULL){
    rep_atom =new int[numCell];
    pdist = new double*[numCell];
    for(int i=0;i<numCell;i++) pdist[i]=new double[3];
  }
  else{
    delete [] rep_atom;
    for(int i=0;i<numCell;i++)delete[] pdist[i];
    delete [] pdist;

    rep_atom = new int[numCell];    
    pdist = new double*[numCell];
    for(int i=0;i<numCell;i++) pdist[i]=new double[3];
  }

  for(int i=0;i<numCell;i++) rep_atom[i]=-1;

  int myrank = get_global_rank();
  
  for(int i=0; i<numPair;i++){
    int rep_ia=-1;
    int rep_ja=-1;
    double dist=-1;
    
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    int inumAtom = listnumAtom[icell];
    int jnumAtom = listnumAtom[jcell];

    if(inumAtom!=0 && jnumAtom!=0){     
      for(int ia=0;ia<inumAtom;ia++){
	for(int ja=0;ja<jnumAtom;ja++){
	  int iatom = atoms[icell][ia];
	  int jatom = atoms[jcell][ja];

	  //double dx = atom->x[iatom][0] - atom->x[jatom][0];
	  //double dy = atom->x[iatom][1] - atom->x[jatom][1];
	  //double dz = atom->x[iatom][2] - atom->x[jatom][2];
	  double dx = gatoms[iatom][0] - gatoms[jatom][0];
	  double dy = gatoms[iatom][1] - gatoms[jatom][1];
	  double dz = gatoms[iatom][2] - gatoms[jatom][2];
	    
	  double l2 = dx*dx + dy*dy + dz*dz;
	  if(ia==0 && ja==0){
	    dist=l2;
	    pdist[jcell][0]=dx;
	    pdist[jcell][1]=dy;
	    pdist[jcell][2]=dz;
	    pdist[icell][0]=-dx;
	    pdist[icell][1]=-dy;
	    pdist[icell][2]=-dz;
	    rep_ia=iatom;
	    rep_ja=jatom;
	  }else if(dist>l2){
	    dist=l2;
	    pdist[jcell][0]=dx;
	    pdist[jcell][1]=dy;
	    pdist[jcell][2]=dz;
	    pdist[icell][0]=-dx;
	    pdist[icell][1]=-dy;
	    pdist[icell][2]=-dz;
	    rep_ia=iatom;
	    rep_ja=jatom;
	  }
	}
      }
    }

    rep_atom[icell]=rep_ia;
    rep_atom[jcell]=rep_ja;

    //cout << icell << " "<<rep_ia <<endl <<atom->x[rep_ia][0] <<" " <<atom->x[rep_ia][1] <<" " <<atom->x[rep_ia][2] <<endl;
    //cout << jcell << " "<<rep_ja <<endl <<atom->x[rep_ja][0] <<" " <<atom->x[rep_ja][1] <<" " <<atom->x[rep_ja][2] <<endl;
    if(myrank==0){
      
      //cout <<"my rank: "<< myrank<<" "<< icell << " "<<rep_ia <<endl <<gatoms[rep_ia][0] <<" " <<gatoms[rep_ia][1] <<" " <<gatoms[rep_ia][2] <<endl;
      //cout << "my rank: "<< myrank<<" " <<jcell << " "<<rep_ja <<endl <<gatoms[rep_ja][0] <<" " <<gatoms[rep_ja][1] <<" " <<gatoms[rep_ja][2] <<endl;
    }
  }

  for(int i=0;i<numPair;i++){
    int icell = pairs[i][0];
    int jcell = pairs[i][1];
    int rep_ia = rep_atom[icell];
    int rep_ja = rep_atom[jcell];
    //cout <<icell << " "<<rep_ia <<" " <<gatoms[rep_ia][0] <<" " <<gatoms[rep_ia][1] <<" " <<gatoms[rep_ia][2] <<endl;
    //cout <<jcell << " "<<rep_ja <<" " <<gatoms[rep_ja][0] <<" " <<gatoms[rep_ja][1] <<" " <<gatoms[rep_ja][2] <<endl;
  }

  
  }*/

void Cells::UpdateLocalRepAtom(){
  UpdateLocal();
  UpdateCell();
  UpdatePair();
  UpdateLocalPair();
  UpdateLocalAtom();

  if(rep_atom==NULL){
    rep_atom =new int[numCell];
    prep_atom =new int[numCell];
    pdist = new double*[numCell];
    sdist = new double*[numCell];
    for(int i=0;i<numCell;i++){
      pdist[i]=new double[3];
      sdist[i]=new double[3];
    }
  }
  else{
    delete [] rep_atom;
    delete [] prep_atom;
    for(int i=0;i<numCell;i++){
      delete[] pdist[i];
      delete[] sdist[i];
    }
    delete [] pdist;
    delete [] sdist;

    rep_atom = new int[numCell];    
    prep_atom = new int[numCell];    
    pdist = new double*[numCell];
    sdist = new double*[numCell];
    for(int i=0;i<numCell;i++){ 
      pdist[i]=new double[3];
      sdist[i]=new double[3];
    }
  }

  for(int i=0;i<numCell;i++) {
    rep_atom[i]=0;
    prep_atom[i]=0;
    pdist[i][0]=0.0;
    pdist[i][1]=0.0;
    pdist[i][2]=0.0;
    sdist[i][0]=0.0;
    sdist[i][1]=0.0;
    sdist[i][2]=0.0;
  }

  //int myrank = get_global_rank();

  if(myrank==0){
    //cout <"The cell pair \n";
    //for(int i=0;i<numPair; i++) cout << lpair[i]<<endl;
  }
  
  for(int i=0; i<numPair;i++){
    int rep_ia=0;
    int rep_ja=0;
    int lrep_ia=0;
    int lrep_ja=0;
    double dist=0.0;
    double dx20=0.0;
    
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    int inumAtom = listnumAtom[icell];
    int jnumAtom = listnumAtom[jcell];
    
    /*for(int ia=0;ia<inumAtom;ia++){
      for(int ja=0;ja<jnumAtom;ja++){
	int iatom = atoms[icell][ia];
	int jatom = atoms[jcell][ja];
	cout <<ia << " " <<iatom <<" "<<gatoms[iatom][0]<<" "<<gatoms[iatom][1]<<" "<<gatoms[iatom][2]<<endl;
	cout <<ja << " " <<jatom <<" "<<gatoms[jatom][0]<<" "<<gatoms[jatom][1]<<" "<<gatoms[jatom][2]<<endl;
      }
      }*/
    
    if(lpair[i]){
      if(inumAtom!=0 && jnumAtom!=0){     
	for(int ia=0;ia<inumAtom;ia++){
	  for(int ja=0;ja<jnumAtom;ja++){
	    int iatom = atoms[icell][ia];
	    int jatom = atoms[jcell][ja];
	    
	    double dx = gatoms[iatom][0] - gatoms[jatom][0];
	    double dy = gatoms[iatom][1] - gatoms[jatom][1];
	    double dz = gatoms[iatom][2] - gatoms[jatom][2];
	    
	    double l2 = dx*dx + dy*dy + dz*dz;
	    double dx2 = dx*dx;
	    if(ia==0 && ja==0){
	      dist=l2;
	      dx20=dx2;//
	      pdist[icell][0]=gatoms[iatom][0];
	      pdist[icell][1]=gatoms[iatom][1];
	      pdist[icell][2]=gatoms[iatom][2];
	      pdist[jcell][0]=gatoms[jatom][0];
	      pdist[jcell][1]=gatoms[jatom][1];
	      pdist[jcell][2]=gatoms[jatom][2];
	      rep_ia=gtags[iatom];
	      rep_ja=gtags[jatom];
	      lrep_ia=iatom;
	      lrep_ja=jatom;
	    }//else if(dist>l2){
	    else if(dx20>dx2){
	      dx20=dx2;
	      dist=l2;
	      pdist[icell][0]=gatoms[iatom][0];
	      pdist[icell][1]=gatoms[iatom][1];
	      pdist[icell][2]=gatoms[iatom][2];
	      pdist[jcell][0]=gatoms[jatom][0];
	      pdist[jcell][1]=gatoms[jatom][1];
	      pdist[jcell][2]=gatoms[jatom][2];
	      rep_ia=gtags[iatom];
	      rep_ja=gtags[jatom];
	      lrep_ia=iatom;
	      lrep_ja=jatom;
	    }
	    /*else if(dx20==dx2){
	      if(dist>l2){
		dx20=dx2;
		dist=l2;
		pdist[icell][0]=gatoms[iatom][0];
		pdist[icell][1]=gatoms[iatom][1];
		pdist[icell][2]=gatoms[iatom][2];
		pdist[jcell][0]=gatoms[jatom][0];
		pdist[jcell][1]=gatoms[jatom][1];
		pdist[jcell][2]=gatoms[jatom][2];
		rep_ia=gtags[iatom];
		rep_ja=gtags[jatom];
		lrep_ia=iatom;
		lrep_ja=jatom;
	      }
	      }*/
	    
	  }
	}
      }
      rep_atom[icell]=rep_ia;
      rep_atom[jcell]=rep_ja;
      prep_atom[icell]=lrep_ia;
      prep_atom[jcell]=lrep_ja;
    }else{
      rep_atom[icell]=rep_ia;
      rep_atom[jcell]=rep_ja;
      prep_atom[icell]=lrep_ia;
      prep_atom[jcell]=lrep_ja;
    }//if
  }

  int *tmp = new int[numCell];
  int *ltmp = new int[numCell];
  double *senddist = new double[3*numCell];
  double *recdist = new double[3*numCell];

  for(int i=0;i<numCell;i++){
    tmp[i]=0;
    ltmp[i]=0;
    senddist[3*i]=pdist[i][0];
    senddist[3*i+1]=pdist[i][1];
    senddist[3*i+2]=pdist[i][2];
    recdist[3*i]=0.0;
    recdist[3*i+1]=0.0;
    recdist[3*i+2]=0.0;
  }
  
  global_isum(rep_atom,tmp,numCell);
  global_isum(prep_atom,ltmp,numCell);
  global_dsum(senddist,recdist,numCell*3);

  for(int i=0;i<numCell;i++){
    rep_atom[i]=tmp[i];
    prep_atom[i]=ltmp[i];
    pdist[i][0]=recdist[3*i];
    pdist[i][1]=recdist[3*i+1];
    pdist[i][2]=recdist[3*i+2];
    sdist[i][0]=recdist[3*i];
    sdist[i][1]=recdist[3*i+1];
    sdist[i][2]=recdist[3*i+2];
  }

  if(myrank==0){
    for(int i=0;i<numPair;i++){
      int icell = pairs[i][0];
      int jcell = pairs[i][1];
      int rep_ia = rep_atom[icell];
      int rep_ja = rep_atom[jcell];
      int lia = prep_atom[icell];
      int lja = prep_atom[jcell];
      cout <<icell << " "<<rep_ia <<" "<<lia <<" " <<pdist[icell][0] <<" " <<pdist[icell][1] <<" " <<pdist[icell][2] <<endl;
      cout <<jcell << " "<<rep_ja <<" " <<lja<<" " <<pdist[jcell][0] <<" " <<pdist[jcell][1] <<" " <<pdist[jcell][2] <<endl;
      fp <<i << " "<<rep_ia <<" " <<lia <<" "<<pdist[icell][0] <<" " <<pdist[icell][1] <<" " <<pdist[icell][2] <<endl;
      fp <<i << " "<<rep_ja <<" " <<lja<<" " <<pdist[jcell][0] <<" " <<pdist[jcell][1] <<" " <<pdist[jcell][2] <<endl;
      
    }
    //fp<<endl;
    
  }
  
  if(seperation){
    dist0 = new double*[numPair]; 
    for(int i=0;i<numPair;i++){
      dist0[i]=new double[3];
    }
    

    for(int i=0;i<numPair;i++){
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    
      for(int j=0;j<3;j++){
	double dx = pdist[icell][j] -pdist[jcell][j];
	dist0[i][j]=fabs(dx);
      
      }
    }

    seperation=false;
  }

  delete [] tmp;
  delete [] ltmp;
  delete [] senddist;
  delete [] recdist;

}

void Cells::UpdateLocalRepStress(){
  UpdateLocal();
  //int nlocal = lammps->atom->nlocal;
  if(nlist==NULL){
    nlist = new int*[numPair]; for(int i=0;i<numPair;i++)nlist[i]=NULL;
    nnum = new int[numPair];
  }
  else{
    delete [] nnum;
    for(int i=0;i<numPair;i++) if(nlist[i]!=NULL)delete[] nlist[i];
    delete [] nlist;

    nlist = new int*[numPair]; for(int i=0;i<numPair;i++)nlist[i]=NULL;
    nnum = new int[numPair];
  }  
  for(int i=0;i<numPair;i++)nnum[i]=0;

  double rcut2 = rcut*rcut;

  double *tmpdist=new double[numCell*3];
  for(int i=0;i<numCell*3;i++)tmpdist[i]=0.0;
  for(int i=0;i<numCell;i++)prep_atom[i]=0;

  for(int i=0; i<numPair;i++){
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    int rep_ia = rep_atom[icell];
    int rep_ja = rep_atom[jcell];

    if(lpair[i]){
      for(int j = 0;j<nlocal;j++){
	if(rep_ia==gtags[j]){
	  prep_atom[icell] = j; 
	  tmpdist[3*icell]=gatoms[j][0];
	  tmpdist[3*icell+1]=gatoms[j][1];
	  tmpdist[3*icell+2]=gatoms[j][2];

	}
	if(rep_ja==gtags[j]){ 
	  prep_atom[jcell] = j; 
	  tmpdist[3*jcell]=gatoms[j][0];
	  tmpdist[3*jcell+1]=gatoms[j][1];
	  tmpdist[3*jcell+2]=gatoms[j][2];
	}
      }
    }
  }

  //cout <<"check1 " <<endl;

  int *ltmp=new int[numCell];
  double *recdist=new double[numCell*3];
  for(int i=0;i<numCell;i++)ltmp[i]=0;
  for(int i=0;i<numCell*3;i++)recdist[i]=0.0;
  global_isum(prep_atom,ltmp,numCell);
  global_dsum(tmpdist,recdist,numCell*3);
  for(int i=0;i<numCell;i++)prep_atom[i]=ltmp[i];
  delete[] ltmp;
  delete[] tmpdist;
  for(int i=0;i<numCell;i++){
    //cout << prep_atom[i] <<" " <<rep_atom[i]<<" " << gtags[prep_atom[i]]<<endl;
  }

  for(int i=0; i<numPair;i++){
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    //if(lpair[i]){
    //cout <<i<<" " << prep_atom[icell] <<" " <<rep_atom[i]<<" " << gtags[prep_atom[i]]<<endl;
    //cout <<i<<" "<< prep_atom[jcell] <<" " <<rep_atom[i]<<" " << gtags[prep_atom[i]]<<endl;
      //}
  }

  for(int i=0; i<numPair;i++){
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    if(lpair[i]){
      //double xi = pdist[icell][0];
      //double yi = pdist[icell][1];
      //double zi = pdist[icell][2];
      int lia= prep_atom[icell];
      
      if(lia >nlocal) lia=0;
      double xi = gatoms[lia][0];
      double yi = gatoms[lia][1];
      double zi = gatoms[lia][2];
      
      vector <int> neigh;
      for(int ia=0;ia<nlocal;ia++){
	double dx = xi-gatoms[ia][0];
	double dy = yi-gatoms[ia][1];
	double dz = zi-gatoms[ia][2];
	double r2 = dx*dx + dy*dy +dz*dz;
	//cout <<r2<<endl;
	if(r2 < rcut2) neigh.push_back(ia);
      }
      
      int nsize = neigh.size();
      nnum[i]=nsize;
      nlist[i]=new int[nsize];
      //cout <<nnum[i]<<endl;
      for(int ni=0;ni<nsize;ni++){ 
	nlist[i][ni] = neigh[ni];
	//cout <<nlist[i][ni]<<" ";
      }
      //cout <<endl;
      
    }
  }


  //cout <<" check1 "<< endl;

  int icompute = lammps->modify->find_compute("voro");
  Compute *voro = lammps->modify->compute[icompute];
  voro->compute_peratom();
  double **dvoro = voro->array_atom;

  icompute = lammps->modify->find_compute("astres");
  Compute *stress = lammps->modify->compute[icompute];
  stress->compute_peratom();
  double **dstress = stress->array_atom;

  //  int nlocal = lammps->atom->nlocal;


  if(repstress==NULL){
    repstress = new double*[numPair];
    avolume =new double[numPair];
  }else{
    for(int i=0;i<numPair;i++) delete [] repstress[i];
    delete [] repstress;
    delete [] avolume;
    repstress = new double*[numPair];
    avolume =new double[numPair];

  }

  //cout <<" check2 "<< endl;
  
  double *sendstress=new double[6*numPair];
  double *recstress=new double[6*numPair];
  double *tmpvolume=new double[numPair];
  
  for(int i=0;i<numPair;i++){
    repstress[i]=new double[6];
    for(int j=0;j<6;j++){
      sendstress[6*i+j]=0.0;
      recstress[6*i+j]=0.0;
      repstress[i][j]=0.0;
    }
    avolume[i]=0.0;
    tmpvolume[i]=0.0;
  }
  
  double lattice = 4.032;
  double svolume = lattice*lattice*lattice/4;
  
  
  for(int i=0; i<numPair;i++){
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    if(lpair[i]){

      int vcount = 0;
      for(int ni=0;ni<nnum[i];ni++){
	int ia=nlist[i][ni];
	//cout << ni <<" " << dvoro[ia][0]<< " " << dvoro[ia][1]<<endl;
	//cout << ni << " "<< dstress[ia][0]<< " " << dstress[ia][1]<<" "<< dstress[ia][2]<<" "<< dstress[ia][3]<<" "<< dstress[ia][4]<<" "<< dstress[ia][5]<<endl;	
	if(svolume*0.8 <dvoro[ia][0] && svolume*1.2>dvoro[ia][0]){
	  avolume[i]+=dvoro[ia][0];
	  vcount++;
	}
	
	for(int si=0;si<6;si++){
	  repstress[i][si]+=dstress[ia][si];
	}
	
      }
      avolume[i]/=(double)vcount;
      
      for(int si=0;si<6;si++){
	repstress[i][si]/=nnum[i];
      }
      double dx = pdist[icell][0]-pdist[jcell][0];
      //cout <<"cell: "<< i <<" dx " << dx << " volume " << avolume[i]<<" stress: " <<(repstress[i][0]+repstress[i][1]+repstress[i][2])/60000.0/avolume[i]<<endl;
    }
  }

  for(int i=0;i<numPair;i++){
    for(int j=0;j<6;j++){
      sendstress[6*i+j]=repstress[i][j];
    }
  }
  
  global_dsum(sendstress,recstress,6*numPair);  
  global_dsum(avolume,tmpvolume,numPair);  
  
  for(int i=0;i<numPair;i++){
    for(int j=0;j<6;j++){
      repstress[i][j]=recstress[6*i+j];
    }
    avolume[i]=tmpvolume[i];
  }

  ofstream frame;
  char frameout[100];
  sprintf(frameout,"tslaw_%d.dat",framenum);
  if(myrank==0)frame.open(frameout);
  framenum++;

  for(int i=0; i<numPair;i++){
    int icell= pairs[i][0];
    int jcell= pairs[i][1];
    int lia = prep_atom[icell];
    int lja = prep_atom[jcell];
    
    int gia = gtags[lia];
    int gja = gtags[lja];

    int pgia = rep_atom[icell];
    int pgja = rep_atom[jcell];
    bool out=true;

    double dx = recdist[3*icell]-recdist[3*jcell];
    double ay = 0.5*(recdist[3*icell+1]+recdist[3*jcell+1]);
    double sx = (pdist[icell][0]-pdist[jcell][0]);
    double d0 = dist0[i][0];
    double sep = fabs(sx)-fabs(d0);

    //if(cbox[icell][3] < crackL ) out=false;
    
    double hydro = (repstress[i][0]+repstress[i][1]+repstress[i][2])/30000.0/avolume[i];
    double pxx = repstress[i][0]/10000.0/avolume[i];
    double pyy = repstress[i][1]/10000.0/avolume[i];
    double pzz = repstress[i][2]/10000.0/avolume[i];

    double pxy = repstress[i][3]/10000.0/avolume[i];
    double pyz = repstress[i][4]/10000.0/avolume[i];
    double pzx = repstress[i][5]/10000.0/avolume[i];

    double vonmise = sqrt(pxx*pxx+pyy*pyy+pzz*pzz - pxx*pyy -pyy*pzz -pzz*pxx +3*(pxy*pxy+pyz*pyz+pzx*pzx));
    double constraint = hydro/vonmise;
      
    if(myrank==0)tsf << i <<" " <<ay << " " << fabs(dx) <<" " <<sep <<" "<< avolume[i]<<" " <<hydro<<" "<< pxx<<" " <<pyy<<" "<<pzz<< " "<< vonmise<<" " << constraint <<" " <<pxy <<" " <<pyz << " "<< pzx << endl;
    if(myrank==0)frame << i <<" " <<ay << " " << fabs(dx) <<" " <<sep <<" "<< avolume[i]<<" " <<hydro<<" "<< pxx<<" " <<pyy<<" "<<pzz<< " "<< vonmise<<" " << constraint<<" " <<pxy <<" " <<pyz << " "<< pzx << endl;
    
  }  

  if(myrank==0)frame.close();

  delete [] recdist;
  delete [] sendstress;
  delete [] recstress;
  delete [] tmpvolume;
  
}


void Cells::UpdateGroup(string s){
  for(int i=0;i<numPair;i++){
    char command1[100]={NULL};
    char command2[100]={NULL};
    char command3[100]={NULL};
    int icell=pairs[i][0];
    int id = rep_atom[icell];
    //cout <<"id: " <<gid <<endl;
    sprintf(command1,"region       \t ba%d sphere %f %f %f 6.0 side in",i,pdist[icell][0],pdist[icell][1],pdist[icell][2]);
    lammps->input->one(command1);    
    sprintf(command2,"group        \t %s%d region ba%d",s.c_str(),i,i);
    lammps->input->one(command2);    
    sprintf(command3,"region       \t ba%d delete",i);
    lammps->input->one(command3);    
  }  
  fp_count++;

  
}

void Cells::UnGroup(string s){
  for(int i=0;i<numPair;i++){
    char command[100]={NULL};
    sprintf(command,"group        \t %s%d delete",s.c_str(),i);
    lammps->input->one(command);    
  }  

}

void Cells::AvgTSLaw(int total, int numavg,string s){
  ifstream fpi;
  ofstream fpo;
  
  myrank=get_global_rank();

  int datanum=9;
  double **data = new double*[numCell];
  for(int i=0;i<numCell;i++){
    data[i]=new double[datanum];
  }
  for(int i=0;i<numCell;i++){
    for(int j=0;j<datanum;j++){
      data[i][j]=0.0;
    }
  }
  
  double *tmp=new double[datanum];


  for(int di=0;di<datanum;di++) tmp[di]=0.0;
  
  int outcount=0;
  int datacount=0;

  int cnt[192][30]={0};
  double value[192][30]={0.0};
  double dist[192][30]={0.0};
  double sepa[30]={0.0};

  for(int i=0;i<30;i++){
    sepa[i] = 20.0/29*i;
  }

  for(int i=0;i<total;i++){
    char filename[100]={NULL};
    sprintf(filename,"%s_%d.dat",s.c_str(),i);
    fpi.open(filename);

    if(fpi.is_open()){
      for(int j=0;j<numCell;j++){
	int temp=0;
	string line;
	getline(fpi,line);
	istringstream(line)>>temp>> tmp[0]>>tmp[1]>>tmp[2]>>tmp[3]>>tmp[4]>>tmp[5]>>tmp[6]>>tmp[7]>>tmp[8];
	for(int di=0;di<4;di++) data[j][di]+=tmp[di];
	for(int di=4;di<datanum;di++) data[j][di]+=2*tmp[di];

	/*for(int si=0;si<29;si++){
	  if( sepa[si]<data[j][1] && data[j][1]<sepa[si+1]){
	    value[j][si]+=data[j][4];//hydro:4 
	    dist[j][si]+=data[j][1];
	    cnt[j][si]++;
	  }
	  }*/
      }
    }    
    fpi.close();

    double standr=25.0;
    if(i%numavg==0){
      char outname[100]={NULL};
      sprintf(outname,"%s_avg_%d.dat",s.c_str(),outcount);
      fpo.open(outname);
      
      for(int j=0;j<numCell;j++){
	if(j>20 && data[j][1]<standr*(numavg+1)) fpo << j << " ";
	cout <<j << " ";
	for(int di=0;di<datanum;di++){
	  data[j][di]/=(double)numavg;
	  cout << data[j][di] <<" ";
	  if(j>20 && data[j][1]<standr*(numavg+1))fpo << data[j][di] <<" ";
	}
	cout <<endl;
	if(j>20&& data[j][1]<standr*(numavg+1))fpo <<endl;
	

      }
      for(int j=0;j<numCell;j++){
	for(int di=0;di<datanum;di++){
	  data[j][di]=0.0;
	}
      }
      fpo.close();      
      outcount++;
      }//if


      /*if(fpi.is_open()){
      for(int j=0;j<numCell;j++){
	int temp=0;
	string line;
	getline(fpi,line);
	istringstream(line)>>temp>> tmp[0]>>tmp[1]>>tmp[2]>>tmp[3]>>tmp[4]>>tmp[5]>>tmp[6]>>tmp[7]>>tmp[8];
	for(int di=0;di<datanum;di++) data[j][di]=tmp[di];
      }
    }    
    fpi.close();
    datacount++;

    char outname[100]={NULL};
    sprintf(outname,"%s_avg_%d.dat",s.c_str(),outcount);
    //if(i%numavg==0){ 
    if(i%2==0){ 
      fpo.open(outname);
      outcount++;
    }

    //if(i%numavg==0){
    if(fpo.is_open()){
      //char outname[100]={NULL};
      //sprintf(outname,"%s_avg_%d.dat",s.c_str(),outcount);
      //fpo.open(outname);
      
      for(int j=20;j<50;j++){
	if( 0.0< data[j][2] && data[j][2]<20) fpo << j << " ";	
	if( 0.0< data[j][2] && data[j][2]<20) cout <<j << " ";
	if( 0.0< data[j][2] && data[j][2]<20){
	  for(int di=0;di<datanum;di++){
	    //data[j][di]/=(double)numavg;
	    cout << data[j][di] <<" ";
	    fpo << data[j][di] <<" ";
	  }
	}
	if( 0.0< data[j][2] && data[j][2]<20) cout <<endl;
	if( 0.0< data[j][2] && data[j][2]<20) fpo <<endl;

      }

      for(int j=0;j<numCell;j++){
	for(int di=0;di<datanum;di++){
	  data[j][di]=0.0;
	}
      }

      //fpo.close();      
      //outcount++;
      }//if
    //if(datacount==(numavg-1)){ 
    //fpo.close();
    //if(datacount==(numavg-1)){ 
    // fpo.close();
    // datacount=0;
    //}*/
    }

  fpi.open("avg.dat");
  double separation[500]={0.0};
  double sevalue[500]={0.0};
  int sepacount[500]={0};
  for(int i=0;i<500;i++)separation[i] = 25.0/449*i;

  if(fpi.is_open()){
    for(int j=0;j<numCell*19;j++){
      int temp=0;
      string line;
      getline(fpi,line);
      istringstream(line)>>temp>> tmp[0]>>tmp[1]>>tmp[2]>>tmp[3]>>tmp[4]>>tmp[5]>>tmp[6]>>tmp[7]>>tmp[8];
      cout<<line;
      for(int i=0;i<499;i++){
	if(separation[i]<tmp[1] && tmp[1]<separation[i+1])  {
	  sevalue[i]+=tmp[5];
	  sepacount[i]++;
	}
      }
    }
  }
  fpi.close();

  fpo.open("hydro.dat");
  for(int i=0;i<500;i++){
    fpo<<separation[i]<<" " <<sevalue[i]/sepacount[i] << " "<< sepacount[i]<<endl;
  }

  fpo.close();

  for(int i=0;i<numCell;i++)delete []data[i];
  delete [] data;
  delete [] tmp;





}

void Cells::AnalyseTSLaw(string s){
  //using namespace std;
  //numPair=16;
  //fp_count=100;

  ifstream *stress = new ifstream[numPair];
  int **numatom = new int*[fp_count];
  double **gstress = new double*[fp_count];

  for(int i=0;i<fp_count;i++) {
    numatom[i]=new int[numPair];
    gstress[i]=new double[6*numPair];
  }
  for(int i=0;i<fp_count;i++) {
    for(int j=0;j<numPair;j++){
      for(int adi=0;adi<6;adi++) gstress[i][6*j+adi]=0.0;
    }
  }
  
  //int myrank=get_global_rank();

  cout <<"Analyse start\n";
  if(myrank==0){
    for(int step=0;step<fp_count;step++){
      for(int i=0;i<numPair;i++){
	ifstream fpi;
	string line;
	char command[100]={NULL};
	sprintf(command,"%s%d_%d.tmp",s.c_str(),i,5000*step);
	
	fpi.open(command);
	
	if(fpi.is_open()){
	  getline(fpi,line);
	  getline(fpi,line);
	  getline(fpi,line);
	  fpi>>numatom[step][i];
	  for(int li = 0;li<6;li++) getline(fpi,line);
	  for(int si=0;si<numatom[step][i];si++){
	    int tmp;
	    double tmp_stres[6]={0.0};
	    //fpi >>tmp;
	    //fpi>> tmp_stres[0]>>tmp_stres[1]>>tmp_stres[2]>>tmp_stres[3]>>tmp_stres[4]>>tmp_stres[5];
	    getline(fpi,line);
	    istringstream(line)>>tmp>> tmp_stres[0]>>tmp_stres[1]>>tmp_stres[2]>>tmp_stres[3]>>tmp_stres[4]>>tmp_stres[5];
	    //cout <<line<<endl;
	    //cout <<tmp <<" " <<tmp_stres[0]<<" " <<tmp_stres[1]<<" " <<tmp_stres[2]<<" " <<tmp_stres[3]<<" " <<tmp_stres[4]<<" " <<tmp_stres[5]<<endl;
	    for(int adi=0;adi<6;adi++){
	      gstress[step][6*i+adi]+=tmp_stres[adi];
	    }

	  }
	  
	  for(int adi=0;adi<6;adi++){
	    gstress[step][6*i+adi]/=(double)numatom[step][i];
	  }
	  
	}else{
	  cout <<"Couldn't open the file "<< command <<endl;
	}
	fpi.close();
      }  
    }
    
    for(int step=0;step<fp_count;step++){
      for(int i=0;i<numPair;i++){
	//cout <<"step :" << step << " group :" << i<< " atom number: " <<numatom[step][i]<<endl;
	//cout <<"stress of group: " << gstress[step][6*i] <<" " <<gstress[step][6*i+1]<<" "<<gstress[step][6*i+2]<<" " <<gstress[step][6*i+3]<<" " <<gstress[step][6*i+4]<<" " <<gstress[step][6*i+5]<<endl;
      }
    }    

  }

  if(myrank==0) fp.close();
  ifstream id;

  if(myrank==0) id.open("id_history.dat");
  double **coordb;
  double **coordt;
  coordb = new double*[fp_count];
  coordt = new double*[fp_count];
  for(int i=0;i<fp_count;i++){
    coordb[i]=new double[3*numPair];
    coordt[i]=new double[3*numPair];
  }
  
  for(int step=0;step<fp_count;step++){
    for(int i=0;i<numPair;i++){
      int tmp1=0, tmp2=0;
      double tmp_dist[3]={0.0};
      string line="                                                    ";
      getline(id,line);
      istringstream(line)>>tmp1>>tmp2>>tmp_dist[0] >> tmp_dist[1]>> tmp_dist[2];
      for(int ix=0;ix<3;ix++) coordb[step][ix] = tmp_dist[ix];
      cout << tmp1<<" " << tmp2<<" " << coordb[step][0]<<" "<<coordb[step][1]<<" "<<coordb[step][2]<<endl;
 
      getline(id,line);
      //cout <<line<<endl;
      istringstream(line)>>tmp1>>tmp2>> tmp_dist[0]>>tmp_dist[1]>>tmp_dist[2];
      for(int ix=0;ix<3;ix++) coordt[step][ix] = tmp_dist[ix];
      cout << tmp1<<" " << tmp2<<" " << coordt[step][0]<<" "<<coordt[step][1]<<" "<<coordt[step][2]<<endl;
    }
    //string line1;
    //getline(id,line1);
    //cout <<endl;
  }
  
  if(myrank==0) id.close();
  
  if(myrank==0)fp.open("TSLaw.dat");
  for(int step=0;step<fp_count;step++){
    for(int i=0;i<numPair;i++){
      double dist = fabs(coordb[step][0] - coordt[step][0]);
      double avg = (gstress[step][6*i]+gstress[step][6*i+1]+gstress[step][6*i+2])/3.0;
      cout <<dist<<" "<<avg<<endl;
      if(myrank==0) fp <<dist<<" "<<avg<<endl;
    }
  }
  if(myrank==0) fp.close();



  

}
