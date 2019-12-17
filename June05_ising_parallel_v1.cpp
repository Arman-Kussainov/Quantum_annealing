#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>	

using namespace std;

int main(int argc, char *argv[]){

  double J=+1,H=0,T,k=1;
  int m=100, n=100;
  int id, ntasks;
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

 // changing value of n to make sure that each process recives the same ammount of data 
  int RowsPerTask = ceil(double(n)/double(ntasks));
  n=RowsPerTask*ntasks;

 // 1D array to collect magnetization data in root process
  double* collect_data = new double [ntasks];

 // 1D spin orientation data recieved by single array
  double* s_recv = new double[m*RowsPerTask];
 // 2D array of spins for indiviudal process for its own chank of data
  double* spin_s = new double[m*RowsPerTask];
  double** s_pins = new double*[m];
 // initialize it
  for (int i = 0; i < m; ++i)
    s_pins[i] = spin_s + RowsPerTask*i;

 // create collective 2D array of spins
  double* s_data = new double[m*n];
  double** s = new double*[m];
  // initialize it
  for (int i = 0; i < m; ++i)
    s[i] = s_data + n*i;
  // write element values
  for(int i=0;i<m;i++){for(int j=0;j<n;j++){s[i][j] = 1.0;}}
  // want to draw a circle
  // double x0=n/2,y0=m/2;
  // for(int i=0;i<m;i++){for(int j=0;j<n;j++){if( pow(double(i)-x0,2)+pow(double(j)-y0,2)<=pow(x0/2,2)  ){s[i][j] = -1.0;}}}

  // create linear array of spins to distribute data between the processes 
  double* s_send = new double[m*n];
  //reassigning 2D data to this 1D array for further distribution
  // COLUMN by COLUMN
   for (int t = 0; t < n; t++){
    for (int q = 0; q < m; q++){
     s_send[t * m + q] = s[q][t];}}

 	MPI_Scatter(s_send,m*RowsPerTask,MPI_DOUBLE,
			s_recv,m*RowsPerTask,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // reconstruct from linear to 2D
   for (int t = 0; t < RowsPerTask; t++){
    for (int q = 0; q < m; q++){
     s_pins[q][t] = s_recv[t * m + q];}}

if(id==0){ofstream file("ising.txt");file.close();
	  ofstream mfile("magnetization.txt");mfile.close();}


///////////////////////////////////////////////////
  int total_steps=5100;  int therm_steps=int(0.2*total_steps);
  double dT=0.01;
 srand ((double(id)+1)*time(NULL));
/////////////////////////////////
  for(double T=2;T<=2.5;T+=dT){
///////////////////////////////
   double M=0;
   for(int c_ount=1;c_ount<=total_steps;c_ount++){

    for (int ii=0;ii<m*RowsPerTask;ii++){    
     int i = rand()%m; int j = rand()%RowsPerTask; 
     int u_p   = i == 0 ? m-1 : i-1;
     int b_ottom  = i == m-1 ? 0 : i+1;
     int l_eft = j == 0 ? RowsPerTask-1 : j-1;
     int r_ight =     j == RowsPerTask-1 ? 0 : j+1;
     double dE=2*J*s_pins[i][j]*(s_pins[u_p][j] + s_pins[b_ottom][j] + s_pins[i][l_eft] + s_pins[i][r_ight]);
     if (exp(-dE/(k*T))>(rand()/(RAND_MAX + 1.0))){s_pins[i][j] = -s_pins[i][j];}
    }
 
    double M_0=0;
    for(int i=0;i<m;i++){for(int j=0;j<RowsPerTask;j++){M_0+=s_pins[i][j];}}
     if (c_ount>=therm_steps){M+=M_0;}
   
 }
 
   cout<<id<<"\t"<<T<<"\t"<<M/double(m*RowsPerTask)/double(0.8*total_steps)<<"\n"; //?? cast double to int

    double M_local=M/double(m*RowsPerTask)/double(0.8*total_steps);
    MPI_Gather(&M_local,1,MPI_DOUBLE,
 	collect_data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // convert from 2D to linear to send back
  for (int t = 0; t < RowsPerTask; t++){
   for (int q = 0; q < m; q++){
    s_recv[t * m + q] = s_pins[q][t];}}

 // send it back to root process
 // This way it will rotate the matrix of final magnetization
    MPI_Gather(s_recv,m*RowsPerTask,MPI_DOUBLE,
 	s_send,m*RowsPerTask,MPI_DOUBLE,0,MPI_COMM_WORLD);

 if(id==0){
  ofstream file("ising.txt",std::ofstream::out | std::ofstream::app);
  double M_total=0;
  for(int i=0;i<ntasks;i++)
   M_total+=collect_data[i];
  file<<T<<"\t"<<M_total/ntasks<<"\n";file.close();
 
  ofstream mfile("magnetization.txt",std::ofstream::out | std::ofstream::app); 
 
   for (int t = 0; t < n; t++){
    for (int q = 0; q < m; q++){
     s[q][t]= s_send[t * m + q];}}

  for (int q = 0; q < m; q++){
   for (int t = 0; t < n; t++){
     mfile<<s[q][t]<<"\t";}mfile<<"\n";}mfile.close();

}
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
