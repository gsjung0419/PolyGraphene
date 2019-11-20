//2013.10.03 G.S. Jung @LAMM
#include "analysis.h"

Analysis::Analysis(){
  numfile=0;
  numline=0;
  numdat=0;
  skipline=0;
  mydata=NULL;
  
}

Analysis::~Analysis(){

}

void Analysis::CloseData(){
  
  for(int i=0; i<numfile;i++) delete  mydata[i];
  delete [] mydata;
  
}

void Analysis::NewData(){
  if(mydata!=NULL){
    cout <<"Error: Check memory allocation!! \n";
    exit(1);
  }
  
  mydata = new double *[numfile];
  for(int i=0; i<numfile;i++){
    mydata[i] = new double[numdat];
  }

  for(int i=0;i<numfile;i++){
    for(int j=0;j<numdat;j++){
      mydata[i][j]=0.0;
    }
  }

}

void Analysis::Average(){
  
  for(int i=0;i<numfile;i++){
    ifstream ifp; 
    char filename[100]={NULL};
    sprintf(filename,"%s_%d.dat",filehead.c_str(),i);
    ifp.open(filename);
    
    //Ignore the first "skipline".. 
    for(int iline=0;iline<skipline;iline++){ 
      string line;      
      getline(ifp,line);
    }
    
    for(int iline=0;iline<numline;iline++){
      string line;
      int temp;
      getline(ifp,line);
      double *tmp = new double[numdat];
      istringstream iss (line);
      iss>>temp;
      for(int idat=0;idat<numdat;idat++){
	iss >>tmp[idat];
      }

      for(int idat=0;idat<numdat;idat++){
	mydata[i][idat]+=tmp[idat];
      }
    }

    for(int idat=0;idat<numdat;idat++){
      mydata[i][idat]/=(double)numline;
    }
    
    ifp.close();
  }

  ofstream ofp;
  ofp.open("average.dat");
  for(int i=0;i<numfile;i++){
    for(int j=0;j<numdat;j++){
      ofp<<mydata[i][j] << " ";
    }
    ofp<<endl;
    
  }

  ofp.close();
}

void Analysis::ThermoLog(){
  ifstream ifp; 
  char filename[100]={NULL};
  sprintf(filename,"%s.dat",filehead.c_str());
  ifp.open(filename);  
  
  ofstream ofp;
  ofp.open("atomicthermo.dat");
  ofp<<"Time1 aVolume2 aPe3 aPress4 aEnthalpy5 atotal6\n";
  
  //Ignore the first "skipline".. 
  for(int iline=0;iline<skipline;iline++){ 
    string line;      
    getline(ifp,line);
  }
  
  int numatom=269;
  for(int iline=0;iline<numline;iline++){
    string line;
    int temp;
    getline(ifp,line);
    double *tmp = new double[numdat];
    istringstream iss (line);
    //iss>>temp;
    for(int idat=0;idat<numdat;idat++){
      iss >>tmp[idat];
    }

    ofp << tmp[0] << " "<< tmp[4]/numatom << " "<<tmp[6]/numatom/23.0609<< " " <<tmp[3]/10000<<" " <<tmp[7]/numatom/23.0609<<" "<<(tmp[5]+tmp[7])/numatom/23.0609<<endl;

  }
  

  
  ifp.close();  
  ofp.close();
  
}
