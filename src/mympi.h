
#ifndef MY_MPI_H
#define MY_MPI_H

#include "precision.h" 

//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END

#include <mpi.h>

void local_isend(int *buf, int count, int dest, int tag);
void local_dsend(RealVar *buf, int count, int dest, int tag);
void local_irecv(int *buf, int count, int source, int tag, MPI_Status *status);
void local_drecv(RealVar *buf, int count, int source, int tag, MPI_Status *status);

//void local_psend(RealVar *buf, int count, MPI_Datatype data_type, int dest, int tag);
//void local_precv(RealVar *buf, int count, MPI_Datatype data_type, int source, int tag, MPI_Status *status);

void global_dsum(RealVar *x, int n);
void global_isum(int *x, int n);
void local_dsum(RealVar *x, int n);
void local_isum(int *x, int n);

int get_local_rank();
int get_local_size();
int get_global_rank();
int get_global_size();
void local_barrier();
void global_barrier();

//extern MPI_Comm mycomm; 

#endif

