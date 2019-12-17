#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>

using namespace std;
const double pi = 3.14159;
int nn_old;

// calculate energy for spin grid with extra edges
// ... deals with its own array IT IS POINTED AT!!

double e_nergy(int *S[], double* cos_theta[], int m, int n, double J){
    double s_um=0;
    for(int i=1; i<m-1; i++){
        for(int j=1; j<n-1; j++){
	s_um=s_um+(cos_theta[S[i][j]][S[i+1][j]]+
			   cos_theta[S[i][j]][S[i-1][j]]+
			   cos_theta[S[i][j]][S[i][j+1]]+
			   cos_theta[S[i][j]][S[i][j-1]]);}}
    return s_um/2*J;}

int main(int argc, char *argv[]){
	int id, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	int qq=2;
	int m=16, n=16; 
	int nn; // will contain the random value for the spin orientation
	
	// the initial array will be expanded to provide each thread with the equal amount of data
	int RowsPerTask = ceil(double(n)/double(ntasks));
	n=RowsPerTask*ntasks;
	
	// now we need to simulate periodic boundary by additional data columns for each process/thread
	// each data piece should get two additional rows
	int N=n+ntasks*2; // two extra column for each data band
	
	// create linear array of spins to distribute data between the processes	
	double* s_send = new double[(m+2)*N];
	
	// 1D array to collect enregy value data in root process
	double* collect_data = new double [ntasks];
 	
	//cout<<setprecision(3);
	double r_atio=0.1;
	double J=1,E1t=0,E2t=0;
	long int total_steps=90000000;
	
	double f=2.71828;
	// Introduce arrays to keep log_gE and H_E data
	int n_odes=m*n*4+1;
	double* log_gE = new double[n_odes];
	double* H_E    = new double[n_odes];
	double* EpN = new double[n_odes];
	for(int i=0;i<n_odes;i++){
		log_gE[i]=0;H_E[i]=0;EpN[i]=double(i-2*m*n);}
		
	// 1D spin orientation data recieved by the single core/thhread
	int* s_recv = new int[(m+2)*(RowsPerTask+2)];

	// 2D array of spins for indiviudal process for its own chank of data
	// PLAYS THE ROLE OF S IN THE SINGLE THREAD PROGRAM
	int* spin_s = new int[(m+2)*(RowsPerTask+2)];
	int** s_pins = new int*[m+2];
	// initialize it
	for(int i = 0; i < m+2; ++i)
		s_pins[i] = spin_s + (RowsPerTask+2)*i;
	
	// ARRAY WITH ORIGINAL DATA
	// create the main, basic 2D array of spins without extra edges
	int* s_data = new int[m*n];
	int** s = new int*[m];
	// initialize it
	for(int i = 0; i < m; i++){
		s[i] = s_data + n*i;}
	// see WikiPedia for Potts model description
	// set element values such that cos(2*pi*s[i][j]/qq) = -1 for qq=2 ==> cos(pi), for qq>2 see by yourself..
	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++){
			s[i][j] = qq-1;}}
			
	//  This array will be holding all original data including the EXTRA boundaries
	int* S_data = new int[(m+2)*(n+2)];
	int** S = new int*[(m+2)];
	for(int i = 0; i < (m+2); ++i){
		S[i] = S_data + (n+2)*i;}
	
	// 2D array of the precomputed cos(theta_j-theta_i) values
	// makes program a little bit faster
	double* cost_heta = new double[qq*qq];
	double** cos_theta = new double*[qq];
	// initialize it
	for(int i = 0; i < qq; ++i){
		cos_theta[i] = cost_heta + qq*i;}
	// precalculating the array
	for(int t = 0; t < qq; t++){
		for (int q = 0; q < qq; q++){
			cos_theta[q][t] = cos(2*pi*double(q-t)/qq);}} // do I need two doubles?
			
	if(id==0){
	// energy of the initial configuration of spins on the grid => E1t
	E1t=e_nergy(S, cos_theta, m+2, n+2, J);}
	
	int c_ount=0;
	while(c_ount<=total_steps){c_ount++;
		for (int ii=0;ii<m*n;ii++){			

			int i = rand()%m; int j = rand()%n; // the indices should fall within the boundaries of s[][]
			int nn =rand()%qq; // nn takes values from 0 to qq-1, because array indexed from 0 ==> cos(2pi*nn/qq)....
			nn_old=s[i][j];s[i][j]=nn;

			int cou_nt=0;
			for(int p_os=1;p_os<=1+(ntasks-1)*(RowsPerTask+2);p_os+=(RowsPerTask+2)){
				for(int j=0+p_os;j<RowsPerTask+p_os;j++){
					for(int i=1;i<(m+2)-1;i++){
						S[i][j] = s[i-1][j-1-2*cou_nt];}}cou_nt++;}
						
			// global top and bottom boundaries
			for(int j=0;j<N;j++){S[0][j]=S[(m+2)-2][j];S[(m+2)-1][j]=S[1][j];}
			// global left and right boundaries
			for(int i=0;i<(m+2);i++){S[i][0]=S[i][N-2];S[i][N-1]=S[i][1];}
			// internal boundaries
			for(int j=0;j<ntasks-1;j++){
				for(int i=0;i<(m+2);i++){
					// innner boundaries between the processes
					S[i][(j+1)*(RowsPerTask+1)+j]=S[i][(j+1)*(RowsPerTask+1)+j+2];
					S[i][(j+1)*(RowsPerTask+1)+j+1]=S[i][(j+1)*(RowsPerTask+1)+j-1];}}
					
			// reassigning 2D data to this 1D array for further distribution
			// COLUMN by COLUMN
			for (int t = 0; t < N; t++){
				for (int q = 0; q < (m+2); q++){
					s_send[t * (m+2) + q] = S[q][t];}}						
			 

			 
		// distribute the data
		// individual process gets its data beyond this point
		MPI_Scatter(s_send,(m+2)*(RowsPerTask+2),MPI_DOUBLE,
		s_recv,(m+2)*(RowsPerTask+2),MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		// reconstruct from linear to 2D in each process
		for (int t = 0; t < RowsPerTask+2; t++){
			for (int q = 0; q < m+2; q++){
				s_pins[q][t] = s_recv[t * (m+2) + q];}}				 
									
		double E_local=e_nergy(s_pins, cos_theta, m+2, RowsPerTask+2, J);
		MPI_Gather(&E_local,1,MPI_DOUBLE,collect_data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);			

		if(id==0){
			
			for(int i=0;i<ntasks;i++)
			E2t+=collect_data[i];
			
			if (exp(log_gE[int(E1t)+2*m*n]-log_gE[int(E2t)+2*m*n])>=rand()/(RAND_MAX + 1.0)){
				E1t=E2t;
				log_gE[int(E2t)+2*m*n]+=log(f);
				H_E[int(E2t)+2*m*n]+=1;
			}else{
				log_gE[int(E1t)+2*m*n]+=log(f);
				H_E[int(E1t)+2*m*n]+=1;
				// no spin configuration update here. Reverse s[i][j] to original value
				s[i][j]=nn_old;}
			}
		}

	if(id==0){
		int max_H=1;
		for(int i=0;i<n_odes;i++){
			if(H_E[i]>=max_H){max_H=H_E[i];}}			
		int min_H=max_H;			
		for(int i=0;i<n_odes;i++){
			if(H_E[i]!=0 && H_E[i]<=min_H){min_H=H_E[i];}}

		r_atio=double(min_H)/double(max_H);  // do I need two doubles?
		//cout<<c_ount<<"\t"<<H_E[0]+H_E[1]<<"\t"<<H_E[4*m*n]<<"\t"<<r_atio<<"\n";		
		if ((H_E[0]+H_E[1]>=0)&&(H_E[4*m*n-1]+H_E[4*m*n]>=0)&&(c_ount>=10000)&&(r_atio>=0.90)){
			cout<<c_ount<<"\t"<<E1t/(m*n)<<"\t"<<min_H<<"\t"<<max_H<<"\t"<<r_atio<<"\t"<<f<<"\n";
			ofstream gfile("log_gE_HE.txt");
			for(int i=0;i<n_odes;i++){
				gfile<<EpN[i]/(m*n)<<"\t"<<log_gE[i]<<"\t"<<H_E[i]<<"\n";}
			gfile.close();
			for(int i=0;i<n_odes;i++){H_E[i]=0;}
			f=sqrt(f);
			c_ount=0;}
		}
				
	}
	
	
	// dynamic arrays clean up
	//************************
	delete[] EpN;
	
	delete[] spin_s;
	delete[] s_pins;
	delete[] s_recv;
	delete[] s_send;

	
	delete[] s_data;
	delete[] s;
	
	delete[] S_data;
	delete[] S;
	
	delete[] cost_heta;
	delete[] cos_theta;

	delete[] collect_data;
	
	delete[] log_gE;
	delete[] H_E;
	
	MPI_Finalize();	
	return 0;
}
