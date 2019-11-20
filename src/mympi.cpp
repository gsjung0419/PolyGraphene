
#include "mympi.h"

extern MPI_Comm myComm;

void global_dsum(RealVar *x, int n)
{
		RealVar *y=new RealVar[n];
		for (int i=0;i<n;i++) {
				y[i]=x[i];
		}
#ifdef SINGLE_PRECISION
		MPI_Allreduce(x,y,n,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
#else
		MPI_Allreduce(x,y,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
		for (int i=0;i<n;i++) {
				x[i]=y[i];
		}
		delete[] y; 
//		*x=*y; 

}

void global_isum(int *x, int n)
{
//		int *y=x;
		int *y=new int[n];
		for (int i=0;i<n;i++) {
				y[i]=x[i];
		}
		MPI_Allreduce(x,y,n,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		for (int i=0;i<n;i++) {
				x[i]=y[i];
		}
		delete[] y; 
//		*x=*y; 
}

void local_dsum(RealVar *x, int n)
{
//		RealVar *y=x;
		RealVar *y=new RealVar[n];
		for (int i=0;i<n;i++) {
//				y[i]=x[i];
				y[i]=0.0;
		}
#ifdef SINGLE_PRECISION
		MPI_Allreduce(x,y,n,MPI_FLOAT,MPI_SUM,myComm);
#else
		MPI_Allreduce(x,y,n,MPI_DOUBLE,MPI_SUM,myComm);
#endif
		for (int i=0;i<n;i++) {
				x[i]=y[i];
		}
		delete[] y; 
//		*x=*y; 
}

void local_isum(int *x, int n)
{
//		int *y=x;
		int *y=new int[n];
		for (int i=0;i<n;i++) {
				y[i]=x[i];
		}
		MPI_Allreduce(x,y,n,MPI_INT,MPI_SUM,myComm);
		for (int i=0;i<n;i++) {
				x[i]=y[i];
		}
		delete[] y; 
//		*x=*y; 
}

void local_isend(int *buf, int count, int dest, int tag)
{
		MPI_Send(buf,count,MPI_INT,dest,tag,myComm);
}

void local_dsend(RealVar *buf, int count, int dest, int tag)
{
#ifdef SINGLE_PRECISION
		MPI_Send(buf,count,MPI_FLOAT,dest,tag,myComm);
#else
		MPI_Send(buf,count,MPI_DOUBLE,dest,tag,myComm);
#endif
}

void local_irecv(int *buf, int count, int source, int tag, MPI_Status *status)
{
		MPI_Recv(buf,count,MPI_INT,source,tag,myComm,status);
}

void local_drecv(RealVar *buf, int count, int source, int tag, MPI_Status *status)
{
#ifdef SINGLE_PRECISION
		MPI_Recv(buf,count,MPI_FLOAT,source,tag,myComm,status);
#else
		MPI_Recv(buf,count,MPI_DOUBLE,source,tag,myComm,status);
#endif
}

/*
void local_psend(RealVar *buf, int count, MPI_Datatype data_type, int dest, int tag)
{
		MPI_Send(buf,count,data_type,dest,tag,myComm);
}

void local_precv(RealVar *buf, int count, MPI_Datatype data_type, int source, int tag, MPI_Status *status)
{
		MPI_Recv(buf,count,data_type,source,tag,myComm,status);
}
*/ 

int get_global_rank()
{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		return rank; 
}

int get_global_size()
{
		int size;
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		return size; 
}


int get_local_rank()
{
		int rank;
		MPI_Comm_rank(myComm,&rank);
		return rank; 
}

int get_local_size()
{
		int size;
		MPI_Comm_size(myComm,&size);
		return size; 
}

void local_barrier()
{
		MPI_Barrier(myComm); 
}

void global_barrier()
{
		MPI_Barrier(MPI_COMM_WORLD); 
}
