#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>	

using namespace std;

int main(int argc, char *argv[]){

  double J=+1,H=0,T,k=1;
  int m=200, n=200;

  int id, ntasks;
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

  int RowsPerTask = ceil(n/ntasks);
  n=RowsPerTask*ntasks;

 // array to collect magnetization data
  double* collect_data = new double [ntasks];

 // create linear recieving arrays
  double* s_recv = new double[m*RowsPerTask];

 // create local 2D array of spins
  double* spin_s = new double[m*n];
  double** s_pins = new double*[m];
 // initialize it
  for (int i = 0; i < m; ++i)
    s_pins[i] = spin_s + RowsPerTask*i;

 // create 2D array of spins
  double* s_data = new double[m*n];
  double** s = new double*[m];
  // initialize it
  for (int i = 0; i < m; ++i)
    s[i] = s_data + n*i;

  // write element values
  for(int i=0;i<m;i++){
   for(int j=0;j<n;j++){
    s[i][j] = 1.0;}}

  // create linear array of spins
  double* s_send = new double[m*n];

  //reassigning to linear array
  for (int q = 0; q < m; q++){
   for (int t = 0; t < n; t++){
   s_send[q * n + t] = s[q][t];}}

 	MPI_Scatter(s_send,m*RowsPerTask,MPI_DOUBLE,
			s_recv,m*RowsPerTask,MPI_DOUBLE,0,MPI_COMM_WORLD);
 
  // collect from linear to 2D
  for (int q = 0; q < m; q++){
   for (int t = 0; t < RowsPerTask; t++){
    s_pins[q][t] = s_recv[q * RowsPerTask + t];}}

if(id==0){ofstream file("ising.txt");file.close();}

//////////////////////////////////////////
  int total_steps=1400;  int therm_steps=int(0.2*total_steps);
  double dT=0.02;
  for(double T=dT;T<=4.5;T+=dT){

   double M=0;
   for(int c_ount=1;c_ount<=total_steps;c_ount++){
    for (int ii=0;ii<m*RowsPerTask;ii++){
    int i = rand()%m; int j = rand()%RowsPerTask; 

    int l_eft   = i == 0 ? m-1 : i-1;
    int r_ight  = i == m-1 ? 0 : i+1;
    int b_ottom = j == 0 ? RowsPerTask-1 : j-1;
    int u_p =     j == RowsPerTask-1 ? 0 : j+1;

    double dE=2*J*s_pins[i][j]*(s_pins[l_eft][j] + s_pins[r_ight][j] + s_pins[i][b_ottom] + s_pins[i][u_p]);
   
    if (exp(-dE/(k*T))>(rand()/(RAND_MAX + 1.0))){s_pins[i][j] = -s_pins[i][j];}
    }
   
    double M_0=0;  for(int i=0;i<m;i++){for(int j=0;j<RowsPerTask;j++){M_0+=s_pins[i][j];}}
 
    if (c_ount>=therm_steps){M+=M_0;}
    }
    cout<<id<<"\t"<<T<<"\t"<<M/double(m*RowsPerTask)/double(0.8*total_steps)<<"\n"; //?? cast double to int

    double M_local=M/double(m*RowsPerTask)/double(0.8*total_steps);
    MPI_Gather(&M_local,1,MPI_DOUBLE,
 	collect_data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

if(id==0){
 ofstream file("ising.txt",std::ofstream::out | std::ofstream::app);
 double M_total=0;
 for(int i=0;i<ntasks;i++)
 M_total+=collect_data[i];
 file<<T<<"\t"<<M_total/ntasks<<"\n";
file.close();}
}
delete[] spin_s;
delete[] s_pins;
delete[] s_recv;

delete[] s_data;
delete[] s;
delete[] s_send;

MPI_Finalize();

  return 0;
}
